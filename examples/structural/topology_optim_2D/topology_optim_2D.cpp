/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2018  Manav Bhatia
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */


// MAST includes
#include "examples/structural/topology_optim_2D/topology_optim_2D.h"
#include "examples/base/input_wrapper.h"
#include "level_set/level_set_discipline.h"
#include "level_set/level_set_system_initialization.h"
#include "level_set/level_set_transient_assembly.h"
#include "level_set/level_set_nonlinear_implicit_assembly.h"
#include "level_set/level_set_reinitialization_transient_assembly.h"
#include "level_set/level_set_volume_output.h"
#include "level_set/level_set_boundary_velocity.h"
#include "elasticity/structural_nonlinear_assembly.h"
#include "elasticity/stress_output_base.h"
#include "elasticity/structural_system_initialization.h"
#include "base/parameter.h"
#include "base/field_function_base.h"
#include "base/mesh_field_function.h"
#include "base/nonlinear_system.h"
#include "base/transient_assembly.h"
#include "solver/first_order_newmark_transient_solver.h"
#include "property_cards/element_property_card_2D.h"

// libMesh includes
#include "libmesh/serial_mesh.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/dof_map.h"
#include "libmesh/exodusII_io.h"


namespace MAST {
    
    class Phi:
    public MAST::FieldFunction<RealVectorX> {
        
    public:
        Phi(Real l1,
            Real l2,
            Real r):
        MAST::FieldFunction<RealVectorX>("Phi"),
        _l1  (l1),
        _l2  (l2),
        _r   (r),
        _pi  (acos(-1.)) {
            
            libmesh_assert_less(r, _l1*.5);
            libmesh_assert_less(r, _l2*.5);
        }
        virtual ~Phi() {}
        virtual void operator()(const libMesh::Point& p,
                                const Real t,
                                RealVectorX& v) const {
            
            libmesh_assert_less_equal(t, 1);
            libmesh_assert_equal_to(v.size(), 1);
            
            
            // circle
            //v(0) = -(pow(p(0)-_l1*.5, 2) + pow(p(1)-_l2*.5, 2) - pow(_r, 2));
            
            // waves
            Real
            nx = 4.,
            ny = 4.,
            c  = 0.5,
            pi = acos(-1.),
            x  = p(0)-.5*_l1,
            y  = p(1)-.5*_l2,
            r  = pow(pow(x,2)+pow(y,2),.5);
            //v(0) = 1.*cos(nx*r*_pi/_l1);
            v(0) = cos(nx*pi*x/_l1)+cos(ny*pi*y/_l2)+c;
            
            // linear
            //v(0) = (p(0)-_l1*0.5)*(-10.);
            //v(0) = (p(0)+p(1)-_l1)*(-10.);
        }
    protected:
        Real
        _l1,
        _l2,
        _r,
        _pi;
    };
    
    class Vel: public MAST::FieldFunction<Real> {
    public:
        Vel(): MAST::FieldFunction<Real>("vel") {}
        
        virtual void operator() (const libMesh::Point& p,
                                 const Real t,
                                 Real& v) const {
            
            // waves
            Real
            nt = 8.,
            th = atan2(p(0)-.15, p(1)-.15);
            v  = sin(nt*th/2.);
            
            // constant
            // v    = 1.;
        }
        
    protected:
        
    };
    
    
    
    class PhiMeshFunction:
    public MAST::FieldFunction<Real> {
    public:
        PhiMeshFunction():
        MAST::FieldFunction<Real>("phi"), _phi(nullptr) { }
        virtual ~PhiMeshFunction(){ if (_phi) delete _phi;}
        
        void init(MAST::SystemInitialization& sys, const libMesh::NumericVector<Real>& sol) {
            if (!_phi) _phi = new MAST::MeshFieldFunction(sys, "phi");
            else _phi->clear();
            _phi->init(sol);
        }
        
        virtual void operator() (const libMesh::Point& p, const Real t, Real& v) const {
            libmesh_assert(_phi);
            RealVectorX v1;
            (*_phi)(p, t, v1);
            v = v1(0);
        }
        
    protected:
        MAST::MeshFieldFunction *_phi;
    };
}


MAST::Examples::TopologyOptimizationLevelSet2D::
TopologyOptimizationLevelSet2D(const libMesh::Parallel::Communicator& comm_in):
MAST::Examples::StructuralExample2D  (comm_in),
MAST::FunctionEvaluation             (comm_in),
_level_set_mesh                      (nullptr),
_level_set_eq_sys                    (nullptr),
_level_set_sys                       (nullptr),
_level_set_sys_init                  (nullptr),
_level_set_discipline                (nullptr),
_level_set_function                  (nullptr),
_level_set_vel                       (nullptr) {
    
}


MAST::Examples::TopologyOptimizationLevelSet2D::~TopologyOptimizationLevelSet2D() {

    if (!_initialized)
        return;
    
    delete _level_set_function;
    delete _level_set_vel;
    delete _level_set_sys_init;
    delete _level_set_discipline;
    delete _level_set_eq_sys;
    delete _level_set_mesh;
    
    for (unsigned int i=0; i<_dv_params.size(); i++)
        delete _dv_params[i].second;
}


void
MAST::Examples::TopologyOptimizationLevelSet2D::initialize_solution() {

    // initialize solution of the structural problem
    MAST::Examples::StructuralExample2D::initialize_solution();
    
    
    // initialize solution of the level set problem
    Real
    length  = (*_input)(_prefix+"length", "length of domain along x-axis", 0.3),
    width   = (*_input)(_prefix+ "width", "length of domain along y-axis", 0.3),
    min_val = std::min(length, width);

    Phi phi(length, width, min_val*0.4);
    _level_set_sys_init->initialize_solution(phi);
}



void
MAST::Examples::TopologyOptimizationLevelSet2D::init(MAST::Examples::GetPotWrapper& input,
                                                     const std::string& prefix) {
    
    // let all other data structures be initialized
    MAST::Examples::StructuralExample2D::init(input, prefix);

    // ask structure to use Mindlin bending operator
    dynamic_cast<MAST::ElementPropertyCard2D&>(*_p_card).set_bending_model(MAST::MINDLIN);
    
    /////////////////////////////////////////////////
    // now initialize the design data.
    /////////////////////////////////////////////////

    // first, initialize the level set functions over the domain
    this->initialize_solution();

    // next, define a new parameter to define design variable for nodal level-set
    // function value
    this->_init_phi_dvs();
    
    // one inequality constraint.
    _n_ineq = 1;
    
    _level_set_vel         = new MAST::LevelSetBoundaryVelocity(2);
    _level_set_function    = new PhiMeshFunction;
}




void
MAST::Examples::TopologyOptimizationLevelSet2D::init_dvar(std::vector<Real>& x,
                                                          std::vector<Real>& xmin,
                                                          std::vector<Real>& xmax) {
    
    // one DV for each element
    x.resize(_n_vars);
    xmin.resize(_n_vars);
    xmax.resize(_n_vars);
    
    std::fill(xmin.begin(), xmin.end(),   -1.);
    std::fill(xmax.begin(), xmax.end(),    1.);
    for (unsigned int i=0; i<_n_vars; i++)
        x[i] = (*_dv_params[i].second)();
}




void
MAST::Examples::TopologyOptimizationLevelSet2D::evaluate(const std::vector<Real>& dvars,
                                                         Real& obj,
                                                         bool eval_obj_grad,
                                                         std::vector<Real>& obj_grad,
                                                         std::vector<Real>& fvals,
                                                         std::vector<bool>& eval_grads,
                                                         std::vector<Real>& grads) {
    
    // copy DVs to level set function
    for (unsigned int i=0; i<_n_vars; i++)
        _level_set_sys->solution->set(_dv_params[i].first, dvars[i]);
    _level_set_sys->solution->close();
    _level_set_function->init(*_level_set_sys_init, *_level_set_sys->solution);

    libMesh::ExodusII_IO(*_level_set_mesh).write_equation_systems("o.exo", *_level_set_eq_sys);

    /**********************************************************************
     * DO NOT zero out the gradient vector, since GCMMA needs it for the  *
     * subproblem solution                                                *
     **********************************************************************/
    MAST::LevelSetNonlinearImplicitAssembly         nonlinear_assembly;
    MAST::LevelSetNonlinearImplicitAssembly         level_set_assembly;
    MAST::StructuralNonlinearAssemblyElemOperations nonlinear_elem_ops;
    nonlinear_elem_ops.set_discipline_and_system(*_discipline, *_structural_sys);
    nonlinear_assembly.set_discipline_and_system(*_discipline, *_structural_sys);
    nonlinear_assembly.set_level_set_function(*_level_set_function);
    level_set_assembly.set_discipline_and_system(*_level_set_discipline, *_level_set_sys_init);
    level_set_assembly.set_level_set_function(*_level_set_function);

    // reinitialize the dof constraints before solution of the linear system
    // FIXME: we should be able to clear the constraint object from the
    // system before it goes out of scope, but libMesh::System does not
    // have a clear method. So, we are going to leave it as is, hoping
    // that libMesh::System will not attempt to use it (most likely, we
    // shoudl be ok).
    _sys->attach_constraint_object(nonlinear_assembly);
    _sys->reinit_constraints();
    
    MAST::LevelSetVolume                            volume(level_set_assembly.get_intersection());
    MAST::StressStrainOutputBase                    stress;
    volume.set_discipline_and_system(*_level_set_discipline, *_level_set_sys_init);
    stress.set_discipline_and_system(*_discipline, *_structural_sys);
    volume.set_participating_elements_to_all();
    stress.set_participating_elements_to_all();
    
    //////////////////////////////////////////////////////////////////////
    // evaluate the objective
    //////////////////////////////////////////////////////////////////////
    level_set_assembly.calculate_output(*_level_set_sys->solution, volume);
    obj       = volume.output_total();
    libMesh::out << "============================\n";

    //////////////////////////////////////////////////////////////////////
    // evaluate the stress constraint
    //////////////////////////////////////////////////////////////////////
    _sys->solve(nonlinear_elem_ops, nonlinear_assembly);
    nonlinear_assembly.calculate_output(*_sys->solution, stress);
    fvals[0]  = stress.output_total();

    //////////////////////////////////////////////////////////////////////
    // evaluate the objective sensitivities, if requested
    //////////////////////////////////////////////////////////////////////
    if (eval_obj_grad)
        _evaluate_volume_sensitivity(volume, level_set_assembly, obj_grad);
    
    //////////////////////////////////////////////////////////////////////
    // check to see if the sensitivity of constraint is requested
    //////////////////////////////////////////////////////////////////////
    bool if_grad_sens = false;
    for (unsigned int i=0; i<eval_grads.size(); i++)
        if_grad_sens = (if_grad_sens || eval_grads[i]);
    
    //////////////////////////////////////////////////////////////////////
    // evaluate the sensitivities for constraints
    //////////////////////////////////////////////////////////////////////
    if (if_grad_sens)
        _evaluate_stress_functional_sensitivity(stress,
                                                nonlinear_elem_ops,
                                                nonlinear_assembly,
                                                eval_grads,
                                                grads);
}




void
MAST::Examples::TopologyOptimizationLevelSet2D::
_evaluate_volume_sensitivity(MAST::LevelSetVolume& volume,
                             MAST::LevelSetNonlinearImplicitAssembly& assembly,
                             std::vector<Real>& obj_grad) {
    
    // iterate over each DV, create a sensitivity vector and calculate the
    // volume sensitivity explicitly
    std::unique_ptr<libMesh::NumericVector<Real>>
    dphi(_level_set_sys->solution->zero_clone().release());
    
    for (unsigned int i=0; i<_n_vars; i++) {
        
        volume.zero();
        dphi->zero();
        dphi->set(_dv_params[i].first, 1.);
        dphi->close();
        _level_set_vel->init(*_level_set_sys_init, *_level_set_sys->solution, *dphi);
        
        assembly.calculate_output_direct_sensitivity(*_level_set_sys->solution,
                                                     *dphi,
                                                     *_dv_params[i].second,
                                                     volume);
        obj_grad[i] = volume.output_sensitivity_total(*_dv_params[i].second);
    }
}



void
MAST::Examples::TopologyOptimizationLevelSet2D::
_evaluate_stress_functional_sensitivity(MAST::StressStrainOutputBase& stress,
                                        MAST::AssemblyElemOperations& elem_ops,
                                        MAST::LevelSetNonlinearImplicitAssembly& assembly,
                                        const std::vector<bool>& eval_grads,
                                        std::vector<Real>& grads) {

    
    _sys->adjoint_solve(elem_ops, stress, assembly, false);

    std::unique_ptr<libMesh::NumericVector<Real>>
    dphi(_level_set_sys->solution->zero_clone().release());
    
    for (unsigned int i=0; i<_n_vars; i++) {
        
        stress.zero();
        dphi->zero();
        dphi->set(_dv_params[i].first, 1.);
        dphi->close();

        // initialize the level set perturbation function to create a velocity
        // field
        _level_set_vel->init(*_level_set_sys_init, *_level_set_sys->solution, *dphi);
        
        assembly.calculate_output_adjoint_sensitivity(*_sys->solution,
                                                      _sys->get_adjoint_solution(),
                                                      *_dv_params[i].second,
                                                      elem_ops,
                                                      stress);
        grads[i] = stress.output_sensitivity_total(*_dv_params[i].second);
    }
}



void
MAST::Examples::TopologyOptimizationLevelSet2D::_init_mesh() {
    
    // first call the parent method
    MAST::Examples::StructuralExample2D::_init_mesh();
    
    _level_set_mesh = new libMesh::SerialMesh(MAST::Examples::StructuralExample2D::comm());
    
    // identify the element type from the input file or from the order
    // of the element
    
    unsigned int
    nx_divs = (*_input)(_prefix+"level_set_nx_divs", "number of elements of level-set mesh along x-axis", 10),
    ny_divs = (*_input)(_prefix+"level_set_ny_divs", "number of elements of level-set mesh along y-axis", 10);
    
    Real
    length  = (*_input)(_prefix+"length", "length of domain along x-axis", 0.3),
    width   = (*_input)(_prefix+ "width", "length of domain along y-axis", 0.3);
    
    std::string
    t = (*_input)(_prefix+"level_set_elem_type", "type of geometric element in the level set mesh", "quad4");
    
    libMesh::ElemType
    e_type = libMesh::Utility::string_to_enum<libMesh::ElemType>(t);
    
    // if high order FE is used, libMesh requires atleast a second order
    // geometric element.
    if (_fetype.order > 1 && e_type == libMesh::QUAD4)
        e_type = libMesh::QUAD9;
    else if (_fetype.order > 1 && e_type == libMesh::TRI3)
        e_type = libMesh::TRI6;
    
    // initialize the mesh with one element
    libMesh::MeshTools::Generation::build_square(*_level_set_mesh,
                                                 nx_divs, ny_divs,
                                                 0, length,
                                                 0, width,
                                                 e_type);
}




void
MAST::Examples::TopologyOptimizationLevelSet2D::_init_eq_sys() {
    
    // first call the parent method
    MAST::Examples::StructuralExample2D::_init_eq_sys();
    _level_set_eq_sys->init();
}




void
MAST::Examples::TopologyOptimizationLevelSet2D::_init_system_and_discipline() {
    
    // first initialize the structural system and discipline
    MAST::Examples::StructuralExample2D::_init_system_and_discipline();

    // FEType to initialize the system
    // get the order and type of element
    std::string
    order_str   = (*_input)(_prefix+ "level_set_fe_order", "order of finite element shape basis functions for level set method",    "first");
    
    libMesh::Order
    o  = libMesh::Utility::string_to_enum<libMesh::Order>(order_str);
    _level_set_fetype = libMesh::FEType(o, libMesh::LAGRANGE);
    

    // now initialize the level set related data structures
    _level_set_eq_sys      = new libMesh::EquationSystems(*_level_set_mesh);
    _level_set_sys         = &(_level_set_eq_sys->add_system<MAST::NonlinearSystem>("level_set"));
    _level_set_sys_init    = new MAST::LevelSetSystemInitialization(*_level_set_sys,
                                                                    _level_set_sys->name(),
                                                                    _level_set_fetype);
    _level_set_discipline  = new MAST::LevelSetDiscipline(*_eq_sys);
}



void
MAST::Examples::TopologyOptimizationLevelSet2D::_init_dirichlet_conditions() {
    
    // constrain only the left and right boundaries
    this->_init_boundary_dirichlet_constraint(1, "right_constraint");
    this->_init_boundary_dirichlet_constraint(3, "left_constraint");
    
    _discipline->init_system_dirichlet_bc(*_sys);
}



void
MAST::Examples::TopologyOptimizationLevelSet2D::_init_loads() {
    
    _init_pressure_load(true, 2);
    //_init_temperature_load();
}



void
MAST::Examples::TopologyOptimizationLevelSet2D::_init_phi_dvs() {

    libmesh_assert(_initialized);
    // this assumes that level set is defined using lagrange shape functions
    libmesh_assert_equal_to(_level_set_fetype.family, libMesh::LAGRANGE);
    
    Real
    tol     = 1.e-6,
    width   = (*_input)(_prefix+ "width", "length of domain along y-axis", 0.3);

    unsigned int
    dof_id  = 0;
    
    Real
    val     = 0.;
    
    // iterate over all the node values
    libMesh::MeshBase::const_node_iterator
    it  = _level_set_mesh->nodes_begin(),
    end = _level_set_mesh->nodes_end();
    
    // maximum number of dvs is the number of nodes on the level set function
    // mesh. We will evaluate the actual number of dvs
    _dv_params.reserve(_level_set_mesh->n_nodes());
    _n_vars = 0;
    
    for ( ; it!=end; it++) {
        
        const libMesh::Node& n = **it;
        
        // only if node is not on the upper edge
        if (std::fabs(n(1)-width) > tol) {
       
            std::ostringstream oss;
            oss << "dv_" << _n_vars;
            dof_id                     = n.dof_number(0, 0, 0);
            val                        = _level_set_sys->solution->el(dof_id);
            
            _dv_params.push_back(std::pair<unsigned int, MAST::Parameter*>());
            _dv_params[_n_vars].first  = dof_id;
            _dv_params[_n_vars].second = new MAST::Parameter(oss.str(), val);
            
            _n_vars++;
        }
    }
}



void
MAST::Examples::TopologyOptimizationLevelSet2D::level_set_solve() {
    
    libmesh_assert(_initialized);
    
    bool
    output      = (*_input)(_prefix+"if_output", "if write output to a file", false),
    propagate   = (*_input)(_prefix+"if_propagate", "if propagate level set, or reinitialize it", true);
    
    
    std::string
    output_name = (*_input)(_prefix+"output_file_root", "prefix of output file names", "output");
    output_name += "_level_set.exo";
    
    // create the nonlinear assembly object
    std::unique_ptr<MAST::TransientAssembly> level_set_assembly;
    if (propagate)
        level_set_assembly.reset(new MAST::TransientAssembly);
    else {
        MAST::LevelSetReinitializationTransientAssembly
        *assembly = new MAST::LevelSetReinitializationTransientAssembly;
        
        libMesh::NumericVector<Real>
        &base_sol = _level_set_sys->add_vector("base_sol");
        base_sol  = *_level_set_sys->solution;
        
        assembly->set_reference_solution(base_sol);
        _level_set_discipline->set_level_set_propagation_mode(false);
        level_set_assembly.reset(assembly);
    }
    MAST::LevelSetTransientAssemblyElemOperations            level_set_elem_ops;
    
    // Transient solver for time integration
    MAST::FirstOrderNewmarkTransientSolver  level_set_solver;
    
    // now solve the system
    level_set_assembly->set_discipline_and_system(*_level_set_discipline,
                                                  *_level_set_sys_init);
    
    // file to write the solution for visualization
    libMesh::ExodusII_IO exodus_writer(*_mesh);
    
    // time solver parameters
    unsigned int
    t_step                         = 0,
    n_steps                        = (*_input)(_prefix+"level_set_n_transient_steps", "number of transient time-steps", 100);
    level_set_solver.dt            = (*_input)(_prefix+"level_set_dt", "time-step size",    1.e-3);
    level_set_solver.beta          = 0.5;
    
    // set the previous state to be same as the current state to account for
    // zero velocity as the initial condition
    level_set_solver.solution(1).zero();
    level_set_solver.solution(1).add(1., level_set_solver.solution());
    level_set_solver.solution(1).close();
    
    
    // loop over time steps
    while (t_step <= n_steps) {
        
        libMesh::out
        << "Time step: " << t_step
        << " :  t = " << _level_set_sys->time
        << std::endl;
        
        // write the time-step
        if (output) {
            
            exodus_writer.write_timestep(output_name,
                                         *_eq_sys,
                                         t_step+1,
                                         _level_set_sys->time);
        }
        
        level_set_solver.solve(*level_set_assembly);
        
        level_set_solver.advance_time_step();
        t_step++;
    }
    
    level_set_assembly->clear_discipline_and_system();
}





//MAST::YoungsModulus::YoungsModulus(const std::string& nm,
//                                   MAST::FieldFunction<Real> &rho,
//                                   const Real base_modulus,
//                                   const Real penalty):
//MAST::FieldFunction<Real>(nm),
//_rho(rho),
//_base_modulus(base_modulus),
//_penalty(penalty)  {
//    
//    _functions.insert(&rho);
//}
//
//
//
//MAST::YoungsModulus::~YoungsModulus() {
//
//}
//
//
//
//void
//MAST::YoungsModulus::operator() (const libMesh::Point& p,
//                                 Real t,
//                                 Real& v) const {
//    
//    _rho(p, t, v);
//    v = _base_modulus*pow(v, _penalty);
//}
//
//
//void
//MAST::YoungsModulus::derivative(   const MAST::FunctionBase& f,
//                                const libMesh::Point& p,
//                                Real t,
//                                Real& v) const {
//    Real
//    r = 0., dr = 0.;
//    
//    _rho(p, t, r);
//    _rho.derivative( f, p, t, dr);
//    v = _base_modulus*pow(r, _penalty-1.)*_penalty*dr;
//}
//
//
//
//
//
//
//
//MAST::TopologyOptimization2D::
//TopologyOptimization2D(const libMesh::Parallel::Communicator& comm):
//MAST::FunctionEvaluation(comm),
//_initialized(false),
//_penalty(0.),
//_volume_fraction(0.),
//_n_divs_x(0),
//_n_divs_y(0),
//_n_elems(0) { }
//
//
//
//
//void
//MAST::TopologyOptimization2D::init(GetPot& infile,
//                                   libMesh::ElemType e_type,
//                                   bool if_vk) {
//    
//    libmesh_assert(!_initialized);
//    
//    // number of elements
//    _n_divs_x    = infile("n_divs_x", 16);
//    _n_divs_y    = infile("n_divs_y", 16);
//    _n_elems     = _n_divs_x*_n_divs_y;
//    if (e_type == libMesh::TRI3)
//        _n_elems    *= 2;
//    
//    
//    // now setup the optimization data
//    _n_vars                = _n_elems; // for thickness variable
//    _n_eq                  = 0;
//    _n_ineq                = 1;        // volume fraction constraint
//    _max_iters             = 1000;
//    
//    
//    
//    // length of domain
//    _length        = infile("length", 0.50);
//    _width         = infile("width",  0.25);
//    
//    // topology optimization parameters
//    _penalty          = infile(        "penalty",  3.);
//    _volume_fraction  = infile("volume_fraction", 0.3);
//    
//    // create the mesh
//    _mesh          = new libMesh::SerialMesh(this->comm());
//    
//    // initialize the mesh with one element
//    libMesh::MeshTools::Generation::build_square(*_mesh,
//                                                 _n_divs_x, _n_divs_y,
//                                                 0, _length,
//                                                 0, _width,
//                                                 e_type);
//    _mesh->write("mesh.exo");
//    // create the equation system
//    _eq_sys    = new  libMesh::EquationSystems(*_mesh);
//    
//    // create the libmesh system
//    _sys       = &(_eq_sys->add_system<MAST::NonlinearSystem>("structural"));
//    _rho_sys   = &(_eq_sys->add_system<libMesh::ExplicitSystem>("density_vars"));
//    
//    // FEType to initialize the system
//    libMesh::FEType fetype     (libMesh::FIRST,    libMesh::LAGRANGE);
//    libMesh::FEType fetype_rho (libMesh::CONSTANT, libMesh::MONOMIAL);
//    
//    
//    // initialize the system to the right set of variables
//    _structural_sys = new MAST::StructuralSystemInitialization(*_sys,
//                                                               _sys->name(),
//                                                               fetype);
//    _discipline     = new MAST::PhysicsDisciplineBase(*_eq_sys);
//
//    // density variable for plotting
//    _rho_sys->add_variable("rho", fetype_rho);
//    
//    // create and add the boundary condition and loads
//    _dirichlet_right  = new MAST::DirichletBoundaryCondition;
//    _dirichlet_left   = new MAST::DirichletBoundaryCondition;
//    
//    _dirichlet_right->init  (1, _structural_sys->vars());
//    _dirichlet_left->init   (3, _structural_sys->vars());
//    
//    _discipline->add_dirichlet_bc(1,  *_dirichlet_right);
//    _discipline->add_dirichlet_bc(3,   *_dirichlet_left);
//    
//    _discipline->init_system_dirichlet_bc(*_sys);
//    
//    // initialize the equation system
//    _eq_sys->init();
//    
//
//    // create the property functions and add them to the
//    _nu              = new MAST::Parameter(  "nu", infile("nu",    0.33));
//    _kappa           = new MAST::Parameter("kappa",infile("kappa",5./6.));
//    _zero            = new MAST::Parameter("zero",                    0.);
//    _press           = new MAST::Parameter(   "p", infile("press", 2.e6));
//    _th              = new MAST::Parameter(  "th", infile("th",   2.e-3));
//    
//    _nu_f            = new MAST::ConstantFieldFunction("nu",          *_nu);
//    _kappa_f         = new MAST::ConstantFieldFunction("kappa",    *_kappa);
//    _hoff_f          = new MAST::ConstantFieldFunction("off",       *_zero);
//    _press_f         = new MAST::ConstantFieldFunction("pressure", *_press);
//    _th_f            = new MAST::ConstantFieldFunction("h",           *_th);
//    
//    // initialize the load
//    _p_load          = new MAST::BoundaryConditionBase(MAST::SURFACE_PRESSURE);
//    _p_load->add(*_press_f);
//    _discipline->add_side_load(2, *_p_load); // pressure on upper surface
//
//    
//    // create the compliance output object
//    _output = new MAST::RealOutputFunction(MAST::STRUCTURAL_COMPLIANCE);
//    _discipline->add_volume_output(0, *_output);
//    
//    
//    // resize the elem vector
//    _elems.resize(_n_elems);
//
//    // element iterators to define property cards
//    libMesh::MeshBase::element_iterator
//    el_it  = _mesh->local_elements_begin(),
//    el_end = _mesh->local_elements_end();
//    
//    unsigned int
//    counter = 0;
//    
//    bool
//    insert_success = false;
//    
//    for (; el_it != el_end; el_it++) {
//        
//        libMesh::Elem* el = *el_it;
//
//        // set the subdomain ID of the element so that the
//        // property card can be associated with this ID
//        el->subdomain_id()      = counter;
//        
//        std::ostringstream oss;
//        oss << "rho_" << counter;
//
//        
//        // create the density parameters for the elements
//        MAST::Parameter*
//        rho    = new MAST::Parameter(oss.str(), 1.);
//        MAST::ConstantFieldFunction*
//        rho_f  = new MAST::ConstantFieldFunction("rho", *rho);
//        MAST::YoungsModulus*
//        E_f    = new MAST::YoungsModulus("E",
//                                         *rho_f,
//                                         70.e9,
//                                         _penalty);
//        
//        
//        // create the material property card
//        MAST::IsotropicMaterialPropertyCard*
//        m_card         = new MAST::IsotropicMaterialPropertyCard;
//        
//        // add the material properties to the card
//        m_card->add(     *E_f);
//        m_card->add(   *_nu_f);
//        m_card->add(   *rho_f);
//        m_card->add(*_kappa_f);
//        
//        // create the element property card
//        MAST::Solid2DSectionElementPropertyCard*
//        p_card         = new MAST::Solid2DSectionElementPropertyCard;
//        
//        // add the section properties to the card
//        p_card->add(*_th_f);
//        p_card->add(*_hoff_f);
//        
//        // tell the section property about the material property
//        p_card->set_material(*m_card);
//        
//        _discipline->set_property_for_subdomain(counter, *p_card);
//        
//        
//        // now add the property and card pointers to the map
//        _elems[counter]         =  el;
//        
//        insert_success =
//        _elem_rho.insert(std::pair<const libMesh::Elem*, MAST::Parameter*>
//                         (el, rho)).second;
//        
//        libmesh_assert(insert_success);
//
//        insert_success =
//        _elem_rho_f.insert(std::pair<const libMesh::Elem*, MAST::ConstantFieldFunction*>
//                           (el, rho_f)).second;
//        
//        libmesh_assert(insert_success);
//
//        
//        insert_success =
//        _elem_E_f.insert(std::pair<const libMesh::Elem*, MAST::YoungsModulus*>
//                         (el, E_f)).second;
//        
//        libmesh_assert(insert_success);
//
//
//        insert_success =
//        _elem_m_card.insert(std::pair<const libMesh::Elem*, MAST::IsotropicMaterialPropertyCard*>
//                            (el, m_card)).second;
//        
//        libmesh_assert(insert_success);
//
//        
//        insert_success =
//        _elem_p_card.insert(std::pair<const libMesh::Elem*, MAST::Solid2DSectionElementPropertyCard*>
//                            (el, p_card)).second;
//        
//        libmesh_assert(insert_success);
//
//        counter++;
//    }
//
//    
//    // create the assembly object
//    _assembly = new MAST::NonlinearImplicitAssembly;
//    _elem_ops = new MAST::StructuralNonlinearAssemblyElemOperations;
//    
//    _assembly->set_discipline_and_system(*_elem_ops,
//                                            *_discipline,
//                                            *_structural_sys);
//    
//    
//    _initialized = true;
//}
//
//
//
//
//MAST::TopologyOptimization2D::~TopologyOptimization2D() {
//    
//    
//    if (_initialized) {
//
//        delete _p_load;
//        delete _dirichlet_right;
//        delete _dirichlet_left;
//        
//        delete _nu_f;
//        delete _kappa_f;
//        delete _hoff_f;
//        delete _press_f;
//        delete _th_f;
//        
//        delete _nu;
//        delete _kappa;
//        delete _zero;
//        delete _press;
//        delete _th;
//        
//        _assembly->clear_discipline_and_system();
//        delete _assembly;
//        delete _elem_ops;
//
//        delete _eq_sys;
//        delete _mesh;
//        
//        delete _discipline;
//        delete _structural_sys;
//        
//        delete _output;
//        
//        // delete the element data
//        {
//            std::map<const libMesh::Elem*, MAST::Parameter*>::iterator
//            it   = _elem_rho.begin(),
//            end  = _elem_rho.end();
//            
//            for (; it != end; it++)
//                delete it->second;
//        }
//
//
//        {
//            std::map<const libMesh::Elem*, MAST::ConstantFieldFunction*>::iterator
//            it   = _elem_rho_f.begin(),
//            end  = _elem_rho_f.end();
//            
//            for (; it != end; it++)
//                delete it->second;
//        }
//
//        
//        {
//            std::map<const libMesh::Elem*, MAST::YoungsModulus*>::iterator
//            it   = _elem_E_f.begin(),
//            end  = _elem_E_f.end();
//            
//            for (; it != end; it++)
//                delete it->second;
//        }
//
//        
//        {
//            std::map<const libMesh::Elem*, MAST::IsotropicMaterialPropertyCard*>::iterator
//            it   = _elem_m_card.begin(),
//            end  = _elem_m_card.end();
//            
//            for (; it != end; it++)
//                delete it->second;
//        }
//
//        
//        {
//            std::map<const libMesh::Elem*, MAST::Solid2DSectionElementPropertyCard*>::iterator
//            it   = _elem_p_card.begin(),
//            end  = _elem_p_card.end();
//            
//            for (; it != end; it++)
//                delete it->second;
//        }
//
//    }
//}
//
//
//
//void
//MAST::TopologyOptimization2D::init_dvar(std::vector<Real>& x,
//                                                          std::vector<Real>& xmin,
//                                                          std::vector<Real>& xmax) {
//    // one DV for each element
//    x.resize(_n_vars);
//    xmin.resize(_n_vars);
//    xmax.resize(_n_vars);
//    
//    std::fill(   x.begin(),    x.end(),   0.3);
//    std::fill(xmin.begin(), xmin.end(), 1.e-5);
//    std::fill(xmax.begin(), xmax.end(),    1.);
//}
//
//
//
//void
//MAST::TopologyOptimization2D::evaluate(const std::vector<Real>& dvars,
//                                                         Real& obj,
//                                                         bool eval_obj_grad,
//                                                         std::vector<Real>& obj_grad,
//                                                         std::vector<Real>& fvals,
//                                                         std::vector<bool>& eval_grads,
//                                                         std::vector<Real>& grads) {
//    
//    
//    libmesh_assert_equal_to(dvars.size(), _n_vars);
//    
//    // set the parameter values equal to the DV value
//    for (unsigned int i=0; i<_n_vars; i++) {
//        const libMesh::Elem* el = _elems[i];
//        MAST::Parameter* par = _elem_rho[el];
//        (*par)() = dvars[i];
//    }
//    
//    // DO NOT zero out the gradient vector, since GCMMA needs it for the
//    // subproblem solution
//    
//    libMesh::Point pt; // dummy point object
//    
//    libMesh::out << "New Eval" << std::endl;
//    
//    // the optimization problem is defined as
//    // min weight, subject to constraints on displacement and stresses
//    
//    
//    
//    //////////////////////////////////////////////////////////////////////
//    // first zero the solution
//    //////////////////////////////////////////////////////////////////////
//    _sys->solution->zero();
//    _output->clear();
//    _sys->solve();
//
//
//    //////////////////////////////////////////////////////////////////////
//    // get the objective and constraints
//    //////////////////////////////////////////////////////////////////////
//    
//    // now get the compliance
//    _assembly->calculate_outputs(*_sys->solution);
//    obj = _output->get_value(); // negative, since J = -K
//    
//    
//    fvals[0] = 0.;
//    Real
//    vol       = 0.,
//    total_vol = 0.;
//    for (unsigned int i=0; i<_n_elems; i++) {
//        vol = _elems[i]->volume();
//        fvals[0]  += dvars[i] * vol; // constraint:  xi vi - V <= 0
//        total_vol += vol;
//    }
//    fvals[0] /= (total_vol * _volume_fraction);
//    fvals[0] -=              1.;
//    
//    
//    
//    //////////////////////////////////////////////////////////////////
//    //   evaluate sensitivity if needed
//    //////////////////////////////////////////////////////////////////
//    
//    // sensitivity of the objective function
//    if (eval_obj_grad) {
//        
//        // we are going to choose to use one parametric sensitivity at a time
//        for (unsigned int i=0; i<_n_elems; i++) {
//            
//            MAST::Parameter* f = _elem_rho[_elems[i]];
//            
//            libMesh::ParameterVector params;
//            params.resize(1);
//            params[0]  = f->ptr();
//            
//            // iterate over each dv and calculate the sensitivity
//            _sys->add_sensitivity_solution(0).zero();
//            
//            // sensitivity analysis
//            _sys->sensitivity_solve(params);
//            _assembly->calculate_output_sensitivity(params,
//                                                    true,
//                                                    *_sys->solution);
//            
//            obj_grad[i] = _output->get_sensitivity(f);
//        }
//    }
//    
//    // now check if the sensitivity of constraint function is requested
//    if (eval_grads[0]) {
//        
//        //////////////////////////////////////////////////////////////////
//        // indices used by GCMMA follow this rule:
//        // grad_k = dfi/dxj  ,  where k = j*NFunc + i
//        //////////////////////////////////////////////////////////////////
//        
//        for (unsigned int i=0; i<_n_elems; i++)
//            grads[i] = _elems[i]->volume()/(_volume_fraction * total_vol);
//    }
//    
//    
//    // write the evaluation output
//    //this->output(0, dvars, obj, fvals, false);
//}
//
//
//
//
//
//
//void
//MAST::TopologyOptimization2D::output(unsigned int iter,
//                                                       const std::vector<Real>& x,
//                                                       Real obj,
//                                                       const std::vector<Real>& fval,
//                                                       bool if_write_to_optim_file) const {
//    
//    libmesh_assert_equal_to(x.size(), _n_vars);
//    
//    // set the desity value in the auxiliary system for output
//    for (unsigned int i=0; i<_n_vars; i++) {
//        const libMesh::Elem* el = _elems[i];
//        _rho_sys->solution->set(el->dof_number(_rho_sys->number(), 0, 0), x[i]);
//    }
//    
//    
//    // write the solution for visualization
//    std::set<std::string> nm;
//    nm.insert(_sys->name());
//    nm.insert(_rho_sys->name());
//    libMesh::ExodusII_IO(*_mesh).write_equation_systems("output.exo",
//                                                        *_eq_sys,
//                                                        &nm);
//    
//    MAST::FunctionEvaluation::output(iter, x, obj, fval, if_write_to_optim_file);
//}
//
//

