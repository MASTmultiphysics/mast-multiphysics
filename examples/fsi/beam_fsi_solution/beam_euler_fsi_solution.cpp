/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2017  Manav Bhatia
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

// C++ includes
#include <iostream>
#include <ostream>


// MAST includes
#include "examples/fsi/beam_fsi_solution/beam_euler_fsi_solution.h"
#include "examples/fluid/meshing/panel_mesh_2D.h"
#include "base/nonlinear_system.h"
#include "fluid/conservative_fluid_system_initialization.h"
#include "fluid/conservative_fluid_discipline.h"
#include "fluid/conservative_fluid_transient_assembly.h"
#include "fluid/pressure_function.h"
#include "fluid/flight_condition.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "base/boundary_condition_base.h"
#include "base/complex_mesh_field_function.h"
#include "elasticity/complex_normal_rotation_mesh_function.h"
#include "elasticity/structural_system_initialization.h"
#include "elasticity/structural_discipline.h"
#include "elasticity/structural_element_base.h"
#include "elasticity/structural_nonlinear_assembly.h"
#include "elasticity/structural_near_null_vector_space.h"
#include "elasticity/normal_rotation_mesh_function.h"
#include "property_cards/solid_1d_section_element_property_card.h"
#include "property_cards/isotropic_material_property_card.h"
#include "boundary_condition/dirichlet_boundary_condition.h"
#include "solver/first_order_newmark_transient_solver.h"
#include "solver/multiphysics_nonlinear_solver.h"


// libMesh includes
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/parameter_vector.h"
#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/nonlinear_solver.h"


extern libMesh::LibMeshInit* __init;


namespace MAST {

    class FSIBoundaryConditionUpdates:
    public MAST::MultiphysicsNonlinearSolverBase::PreResidualUpdate {
        
    public:
        FSIBoundaryConditionUpdates(MAST::SystemInitialization&             structural_sys,
                                    MAST::SystemInitialization&             fluid_sys,
                                    MAST::MeshFieldFunction&                vel,
                                    MAST::MeshFieldFunction&                displ,
                                    MAST::PressureFunction&                 press):
        _structural_sys   (structural_sys),
        _fluid_sys        (fluid_sys),
        _vel              (vel),
        _displ            (displ),
        _press            (press)
        { }
        
        virtual ~FSIBoundaryConditionUpdates() { }
        
        virtual void
        update_at_solution(std::vector<libMesh::NumericVector<Real>*>&  sol_vecs) {

            libmesh_error(); // setup velocity update
            
            // make sure that the solutions are appropriately sized
            libMesh::NumericVector<Real>
            &fluid_sol      = *sol_vecs[0],
            &structural_sol = *sol_vecs[1];
            
            libmesh_assert_equal_to(fluid_sol.size(),
                                    _fluid_sys.system().n_dofs());
            libmesh_assert_equal_to(fluid_sol.local_size(),
                                    _fluid_sys.system().n_local_dofs());

            libmesh_assert_equal_to(structural_sol.size(),
                                    _structural_sys.system().n_dofs());
            libmesh_assert_equal_to(structural_sol.local_size(),
                                    _structural_sys.system().n_local_dofs());

            _displ.init      (structural_sol);
            _press.init      (fluid_sol);
        }

        
        virtual void
        update_at_perturbed_solution(std::vector<libMesh::NumericVector<Real>*>&   sol_vecs,
                                     std::vector<libMesh::NumericVector<Real>*>&  dsol_vecs) {
            
            libmesh_error(); // setup velocity update

            // make sure that the solutions are appropriately sized
            libMesh::NumericVector<Real>
            &fluid_sol           = * sol_vecs[0],
            &fluid_sol_sens      = *dsol_vecs[0],
            &structural_sol      = * sol_vecs[1],
            &structural_sol_sens = *dsol_vecs[1];
            
            libmesh_assert_equal_to(fluid_sol.size(),
                                    _fluid_sys.system().n_dofs());
            libmesh_assert_equal_to(fluid_sol.local_size(),
                                    _fluid_sys.system().n_local_dofs());
            libmesh_assert_equal_to(fluid_sol_sens.size(),
                                    _fluid_sys.system().n_dofs());
            libmesh_assert_equal_to(fluid_sol_sens.local_size(),
                                    _fluid_sys.system().n_local_dofs());
            
            libmesh_assert_equal_to(structural_sol.size(),
                                    _structural_sys.system().n_dofs());
            libmesh_assert_equal_to(structural_sol.local_size(),
                                    _structural_sys.system().n_local_dofs());
            libmesh_assert_equal_to(structural_sol_sens.size(),
                                    _structural_sys.system().n_dofs());
            libmesh_assert_equal_to(structural_sol_sens.local_size(),
                                    _structural_sys.system().n_local_dofs());
            
            //_displ.init(structural_sol);
            //_normal_rot.init(structural_sol);
            _press.init(fluid_sol, &fluid_sol_sens);
        }

        
    protected:
        
        MAST::SystemInitialization  &           _structural_sys;
        MAST::SystemInitialization  &           _fluid_sys;
        MAST::MeshFieldFunction     &           _vel;
        MAST::MeshFieldFunction     &           _displ;
        MAST::PressureFunction      &           _press;
    };
}


MAST::BeamEulerFSIAnalysis::BeamEulerFSIAnalysis():
_structural_mesh                    (nullptr),
_fluid_mesh                         (nullptr),
_structural_eq_sys                  (nullptr),
_fluid_eq_sys                       (nullptr),
_structural_sys                     (nullptr),
_fluid_sys                          (nullptr),
_structural_sys_init                (nullptr),
_structural_discipline              (nullptr),
_fluid_sys_init                     (nullptr),
_fluid_discipline                   (nullptr),
_flight_cond                        (nullptr),
_far_field                          (nullptr),
_symm_wall                          (nullptr),
_slip_wall                          (nullptr),
_pressure                           (nullptr),
_displ                              (nullptr),
_normal_rot                         (nullptr),
_pressure_function                  (nullptr),
_velocity                           (nullptr),
_velocity_f                         (nullptr),
_length                             (0.),
_thy                                (nullptr),
_thz                                (nullptr),
_rho                                (nullptr),
_E                                  (nullptr),
_nu                                 (nullptr),
_zero                               (nullptr),
_mach                               (nullptr),
_rho_air                            (nullptr),
_gamma_air                          (nullptr),
_thy_f                              (nullptr),
_thz_f                              (nullptr),
_rho_f                              (nullptr),
_E_f                                (nullptr),
_nu_f                               (nullptr),
_hyoff_f                            (nullptr),
_hzoff_f                            (nullptr),
_mach_f                             (nullptr),
_rho_air_f                          (nullptr),
_gamma_air_f                        (nullptr),
_m_card                             (nullptr),
_p_card                             (nullptr),
_dirichlet_left                     (nullptr),
_dirichlet_right                    (nullptr),
_bc_updates                         (nullptr) {
    
    //////////////////////////////////////////////////////////////////////
    //    SETUP THE FLUID DATA
    //////////////////////////////////////////////////////////////////////
    
    // initialize the libMesh object
    _fluid_mesh              = new libMesh::ParallelMesh(__init->comm());
    _fluid_eq_sys            = new libMesh::EquationSystems(*_fluid_mesh);
    
    
    // add the system to be used for analysis
    _fluid_sys = &(_fluid_eq_sys->add_system<MAST::NonlinearSystem>("fluid"));
    _fluid_sys->set_init_B_matrix();
    
    
    // initialize the flow conditions
    GetPot infile("input.in");
    
    
    const unsigned int
    dim                 = 2,
    nx_divs             = 3,
    ny_divs             = 1,
    panel_bc_id         = 10,
    symmetry_bc_id      = 11;
    
    libMesh::ElemType
    elem_type           =
    libMesh::Utility::string_to_enum<libMesh::ElemType>(infile("elem_type", "QUAD4"));
    
    libMesh::FEFamily
    fe_type             =
    libMesh::Utility::string_to_enum<libMesh::FEFamily>(infile("fe_family", "LAGRANGE"));
    
    libMesh::Order
    fe_order            =
    libMesh::Utility::string_to_enum<libMesh::Order>(infile("fe_order", "FIRST"));
    
    std::vector<Real>
    x_div_loc        (nx_divs+1),
    x_relative_dx    (nx_divs+1),
    y_div_loc        (ny_divs+1),
    y_relative_dx    (ny_divs+1);
    
    std::vector<unsigned int>
    x_divs           (nx_divs),
    y_divs           (ny_divs);
    
    std::auto_ptr<MAST::MeshInitializer::CoordinateDivisions>
    x_coord_divs    (new MAST::MeshInitializer::CoordinateDivisions),
    y_coord_divs    (new MAST::MeshInitializer::CoordinateDivisions);
    
    std::vector<MAST::MeshInitializer::CoordinateDivisions*>
    divs(dim);
    
    
    // now read in the values: x-coord
    for (unsigned int i_div=0; i_div<nx_divs+1; i_div++) {
        
        x_div_loc[i_div]        = infile("x_div_loc",   0., i_div);
        x_relative_dx[i_div]    = infile( "x_rel_dx",   0., i_div);
        
        if (i_div < nx_divs) //  this is only till nx_divs
            x_divs[i_div]       = infile( "x_div_nelem", 0, i_div);
    }
    
    divs[0] = x_coord_divs.get();
    x_coord_divs->init(nx_divs, x_div_loc, x_relative_dx, x_divs);
    
    
    // now read in the values: y-coord
    for (unsigned int i_div=0; i_div<ny_divs+1; i_div++) {
        
        y_div_loc[i_div]     = infile("y_div_loc", 0., i_div);
        y_relative_dx[i_div] = infile( "y_rel_dx", 0., i_div);
        
        if (i_div < ny_divs) //  this is only till ny_divs
            y_divs[i_div]    = infile( "y_div_nelem",  0, i_div);
    }
    
    divs[1] = y_coord_divs.get();
    y_coord_divs->init(ny_divs, y_div_loc, y_relative_dx, y_divs);
    
    
    
    
    // initialize the mesh
    MAST::PanelMesh2D().init(0.,               // t/c
                             false,            // if cos bump
                             0,                // n max bumps
                             panel_bc_id,
                             symmetry_bc_id,
                             divs,
                             *_fluid_mesh,
                             elem_type);
    
    _fluid_discipline   = new MAST::ConservativeFluidDiscipline(*_fluid_eq_sys);
    _fluid_sys_init     = new MAST::ConservativeFluidSystemInitialization(*_fluid_sys,
                                                                          _fluid_sys->name(),
                                                                          libMesh::FEType(fe_order, fe_type),
                                                                          dim);
    
    
    // initialize the equation system for analysis
    _fluid_eq_sys->init();
    
    // print the information
    _fluid_eq_sys->print_info();
    
    // create the oundary conditions for slip-wall and far-field
    _far_field     = new MAST::BoundaryConditionBase(MAST::FAR_FIELD),
    _symm_wall     = new MAST::BoundaryConditionBase(MAST::SYMMETRY_WALL);
    _slip_wall     = new MAST::BoundaryConditionBase(MAST::SLIP_WALL);
    
    _flight_cond    =  new MAST::FlightCondition;
    for (unsigned int i=0; i<3; i++) {
        
        _flight_cond->body_roll_axis(i)     = infile(    "body_roll_axis", 0., i);
        _flight_cond->body_pitch_axis(i)    = infile(   "body_pitch_axis", 0., i);
        _flight_cond->body_yaw_axis(i)      = infile(     "body_yaw_axis", 0., i);
        _flight_cond->body_euler_angles(i)  = infile( "body_euler_angles", 0., i);
        _flight_cond->body_angular_rates(i) = infile("body_angular_rates", 0., i);
    }
    
    _flight_cond->ref_chord       = infile("ref_c",    1.);
    _flight_cond->altitude        = infile( "alt",     0.);
    _flight_cond->mach            = infile("mach",     .5);
    _flight_cond->gas_property.cp = infile(  "cp",  1003.);
    _flight_cond->gas_property.cv = infile(  "cv",   716.);
    _flight_cond->gas_property.T  = infile("temp",   300.);
    _flight_cond->gas_property.rho= infile( "rho",   1.05);
    
    _flight_cond->init();
    
    // tell the discipline about the fluid values
    _fluid_discipline->set_flight_condition(*_flight_cond);
    
    // define parameters
    _velocity          = new MAST::Parameter("velocity",  _flight_cond->velocity_magnitude);
    
    
    // now define the constant field functions based on this
    _velocity_f        = new MAST::ConstantFieldFunction("velocity", *_velocity);
    
    // tell the physics about boundary conditions
    _fluid_discipline->add_side_load(    panel_bc_id, *_slip_wall);
    _fluid_discipline->add_side_load( symmetry_bc_id, *_symm_wall);
    // all boundaries except the bottom are far-field
    for (unsigned int i=1; i<=3; i++)
        _fluid_discipline->add_side_load(              i, *_far_field);
    
    _pressure_function =
    new MAST::PressureFunction(*_fluid_sys_init, *_flight_cond);
    
        
    //////////////////////////////////////////////////////////////////////
    //    SETUP THE STRUCTURAL DATA
    //////////////////////////////////////////////////////////////////////
    
    x_div_loc.resize     (2);
    x_relative_dx.resize (2);
    x_divs.resize        (1);
    divs.resize          (1);
    
    x_coord_divs.reset   (new MeshInitializer::CoordinateDivisions);
    
    
    // now read in the values: x-coord
    for (unsigned int i_div=0; i_div<2; i_div++) {
        
        x_div_loc[i_div]        = infile("x_div_loc",   0., i_div+1);
        x_relative_dx[i_div]    = infile( "x_rel_dx",   0., i_div+1);
    }
    x_divs[0]       = infile( "x_div_nelem", 0, 1);
    
    divs[0] = x_coord_divs.get();
    x_coord_divs->init(1, x_div_loc, x_relative_dx, x_divs);
    
    
    // setup length for use in setup of flutter solver
    _length = x_div_loc[1]-x_div_loc[0];
    
    // create the mesh
    _structural_mesh       = new libMesh::SerialMesh(__init->comm());
    
    
    MeshInitializer().init(divs, *_structural_mesh, libMesh::EDGE2);
    
    // create the equation system
    _structural_eq_sys    = new  libMesh::EquationSystems(*_structural_mesh);
    
    // create the libmesh system
    _structural_sys       = &(_structural_eq_sys->add_system<MAST::NonlinearSystem>("structural"));
    _structural_sys->set_eigenproblem_type(libMesh::GHEP);
    
    // FEType to initialize the system
    libMesh::FEType fetype (libMesh::FIRST, libMesh::LAGRANGE);
    
    // initialize the system to the right set of variables
    _structural_sys_init  = new MAST::StructuralSystemInitialization(*_structural_sys,
                                                                     _structural_sys->name(),
                                                                     fetype);
    _structural_discipline = new MAST::StructuralDiscipline(*_structural_eq_sys);
    
    
    // create and add the boundary condition and loads
    _dirichlet_left = new MAST::DirichletBoundaryCondition;
    _dirichlet_right= new MAST::DirichletBoundaryCondition;
    std::vector<unsigned int> constrained_vars(4);
    constrained_vars[0] = 0;  // u
    constrained_vars[1] = 1;  // v
    constrained_vars[2] = 2;  // w
    constrained_vars[3] = 3;  // tx
    _dirichlet_left->init (0, constrained_vars);
    _dirichlet_right->init(1, constrained_vars);
    _structural_discipline->add_dirichlet_bc(0, *_dirichlet_left);
    _structural_discipline->add_dirichlet_bc(1, *_dirichlet_right);
    _structural_discipline->init_system_dirichlet_bc(*_structural_sys);
    
    // initialize the equation system
    _structural_eq_sys->init();
    
    // initialize the motion object
    _vel          = new MAST::MeshFieldFunction(*_structural_sys_init,
                                                "velocity");
    _displ        = new MAST::MeshFieldFunction(*_structural_sys_init,
                                                "displacement");
    _normal_rot   = new MAST::NormalRotationMeshFunction("frequency_domain_normal_rotation",
                                                         *_displ);
    _slip_wall->add(*_vel);
    _slip_wall->add(*_normal_rot);
    
    
    _structural_sys->eigen_solver->set_position_of_spectrum(libMesh::LARGEST_MAGNITUDE);
    _structural_sys->set_exchange_A_and_B(true);
    _structural_sys->set_n_requested_eigenvalues(infile("n_modes", 10));
    
    // create the property functions and add them to the
    
    _thy             = new MAST::Parameter("thy",  0.0015);
    _thz             = new MAST::Parameter("thz",    1.00);
    _rho             = new MAST::Parameter("rho",   2.7e3);
    _E               = new MAST::Parameter("E",     72.e9);
    _nu              = new MAST::Parameter("nu",     0.33);
    _zero            = new MAST::Parameter("zero",     0.);
    _mach            = new MAST::Parameter("mach",     3.);
    _rho_air         = new MAST::Parameter("rho" ,   1.05);
    _gamma_air       = new MAST::Parameter("gamma",   1.4);
    
    
    
    // prepare the vector of parameters with respect to which the sensitivity
    // needs to be benchmarked
    _params_for_sensitivity.push_back(_E);
    _params_for_sensitivity.push_back(_nu);
    _params_for_sensitivity.push_back(_thy);
    _params_for_sensitivity.push_back(_thz);
    
    
    
    _thy_f           = new MAST::ConstantFieldFunction("hy",          *_thy);
    _thz_f           = new MAST::ConstantFieldFunction("hz",          *_thz);
    _rho_f           = new MAST::ConstantFieldFunction("rho",         *_rho);
    _E_f             = new MAST::ConstantFieldFunction("E",             *_E);
    _nu_f            = new MAST::ConstantFieldFunction("nu",           *_nu);
    _hyoff_f         = new MAST::ConstantFieldFunction("hy_off",     *_zero);
    _hzoff_f         = new MAST::ConstantFieldFunction("hz_off",     *_zero);
    _mach_f          = new MAST::ConstantFieldFunction("mach",       *_mach);
    _rho_air_f       = new MAST::ConstantFieldFunction("rho",     *_rho_air);
    _gamma_air_f     = new MAST::ConstantFieldFunction("gamma", *_gamma_air);
    
    // create the material property card
    _m_card          = new MAST::IsotropicMaterialPropertyCard;
    
    // add the material properties to the card
    _m_card->add(*_rho_f);
    _m_card->add(*_E_f);
    _m_card->add(*_nu_f);
    
    // create the element property card
    _p_card          = new MAST::Solid1DSectionElementPropertyCard;
    
    // tell the card about the orientation
    libMesh::Point orientation;
    orientation(1) = 1.;
    _p_card->y_vector() = orientation;
    
    // add the section properties to the card
    _p_card->add(*_thy_f);
    _p_card->add(*_thz_f);
    _p_card->add(*_hyoff_f);
    _p_card->add(*_hzoff_f);
    
    // tell the section property about the material property
    _p_card->set_material(*_m_card);
    
    _p_card->init();
    
    _structural_discipline->set_property_for_subdomain(0, *_p_card);
    
    // pressure boundary condition for the beam
    _pressure    =  new MAST::BoundaryConditionBase(MAST::SURFACE_PRESSURE);
    _pressure->add(*_pressure_function);
    _structural_discipline->add_volume_load(0, *_pressure);

    
}






MAST::BeamEulerFSIAnalysis::~BeamEulerFSIAnalysis() {
    
    delete _fluid_eq_sys;
    delete _fluid_mesh;
    
    delete _fluid_discipline;
    delete _fluid_sys_init;
    
    delete _far_field;
    delete _symm_wall;
    delete _slip_wall;
    
    delete _flight_cond;
    
    delete _velocity;
    
    delete _velocity_f;
    
    delete _pressure_function;    
    
    delete _displ;
    delete _normal_rot;
    delete _pressure;
    
    delete _structural_eq_sys;
    delete _structural_mesh;
    
    
    delete _structural_discipline;
    delete _structural_sys_init;
    
    delete _m_card;
    delete _p_card;
    
    delete _dirichlet_left;
    delete _dirichlet_right;
    
    delete _thy_f;
    delete _thz_f;
    delete _rho_f;
    delete _E_f;
    delete _nu_f;
    delete _hyoff_f;
    delete _hzoff_f;
    delete _mach_f;
    delete _rho_air_f;
    delete _gamma_air_f;
    
    
    delete _thy;
    delete _thz;
    delete _rho;
    delete _E;
    delete _nu;
    delete _zero;
    delete _mach;
    delete _rho_air;
    delete _gamma_air;
    
    delete _bc_updates;
}





MAST::Parameter*
MAST::BeamEulerFSIAnalysis::get_parameter(const std::string &nm) {
    
    MAST::Parameter *rval = nullptr;
    
    // look through the vector of parameters to see if the name is available
    std::vector<MAST::Parameter*>::iterator
    it   =  _params_for_sensitivity.begin(),
    end  =  _params_for_sensitivity.end();
    
    bool
    found = false;
    
    for ( ; it != end; it++) {
        
        if (nm == (*it)->name()) {
            rval    = *it;
            found   = true;
        }
    }
    
    // if the param was not found, then print the message
    if (!found) {
        libMesh::out
        << std::endl
        << "Parameter not found by name: " << nm << std::endl
        << "Valid names are: "
        << std::endl;
        for (it = _params_for_sensitivity.begin(); it != end; it++)
            libMesh::out << "   " << (*it)->name() << std::endl;
        libMesh::out << std::endl;
    }
    
    return rval;
}





Real
MAST::BeamEulerFSIAnalysis::solve(bool if_write_output,
                                  const Real tol,
                                  const unsigned int max_bisection_iters) {
    
    
    /////////////////////////////////////////////////////////////////
    //  INITIALIZE FLUID SOLUTION
    /////////////////////////////////////////////////////////////////
    // the modal and flutter problems are solved on rank 0, while
    // the fluid solution is setup on the global communicator
    
    // initialize the solution
    RealVectorX s = RealVectorX::Zero(4);
    s(0) = _flight_cond->rho();
    s(1) = _flight_cond->rho_u1();
    s(2) = _flight_cond->rho_u2();
    s(3) = _flight_cond->rho_e();
    _fluid_sys_init->initialize_solution(s);
    
    // create the nonlinear assembly object
    MAST::ConservativeFluidTransientAssembly        fluid_assembly;
    MAST::FirstOrderNewmarkTransientSolver          transient_solver;
    
    // now setup the assembly object
    fluid_assembly.attach_discipline_and_system(*_fluid_discipline,
                                                transient_solver,
                                                *_fluid_sys_init);

    // time solver parameters
    unsigned int
    t_step            = 0,
    n_iters_change_dt = 4,
    iter_count_dt     = 0;
    
    Real
    tval       = 0.,
    vel_0      = 0.,
    vel_1      = 1.e+12,
    p          = 0.5,
    factor     = 0.,
    min_factor = 1.5;
    
    transient_solver.dt            = 1.e-4;
    transient_solver.beta          = 1.0;
    
    // set the previous state to be same as the current state to account for
    // zero velocity as the initial condition
    transient_solver.solution(1).zero();
    transient_solver.solution(1).add(1., transient_solver.solution());
    transient_solver.solution(1).close();
    
    
    MAST::StructuralNonlinearAssembly   structural_assembly;
    
    // create the nonlinear assembly object
    structural_assembly.attach_discipline_and_system(*_structural_discipline,
                                                     *_structural_sys_init);
    MAST::StructuralNearNullVectorSpace nsp;
    _structural_sys->nonlinear_solver->nearnullspace_object = &nsp;
    
    ///////////////////////////////////////////////////////////////////
    // FSI SOLUTION
    ///////////////////////////////////////////////////////////////////
    
    MAST::FSIBoundaryConditionUpdates bc_updates(*_structural_sys_init,
                                                 *_fluid_sys_init,
                                                 *_vel,
                                                 *_displ,
                                                 *_pressure_function);
    
    MAST::MultiphysicsNonlinearSolverBase fsi_solver(__init->comm(), "fsi", 2);

    
    fsi_solver.set_pre_residual_update_object(bc_updates);
    fsi_solver.set_system_assembly(0,      fluid_assembly);
    fsi_solver.set_system_assembly(1, structural_assembly);
    
    fsi_solver.solve();
    
    fluid_assembly.clear_discipline_and_system();
    structural_assembly.clear_discipline_and_system();
    
    return 0.;
}





Real
MAST::BeamEulerFSIAnalysis::sensitivity_solve(MAST::Parameter& p) {
    
    libmesh_error();
    
    return 0.;
}

