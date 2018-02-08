///*
// * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
// * Copyright (C) 2013-2018  Manav Bhatia
// *
// * This library is free software; you can redistribute it and/or
// * modify it under the terms of the GNU Lesser General Public
// * License as published by the Free Software Foundation; either
// * version 2.1 of the License, or (at your option) any later version.
// *
// * This library is distributed in the hope that it will be useful,
// * but WITHOUT ANY WARRANTY; without even the implied warranty of
// * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// * Lesser General Public License for more details.
// *
// * You should have received a copy of the GNU Lesser General Public
// * License along with this library; if not, write to the Free Software
// * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
// */
//
//// C++ includes
//#include <iostream>
//#include <numeric>
//
//
//// MAST includes
//#include "examples/structural/plate_bending_level_set/plate_bending_level_set.h"
//#include "elasticity/structural_system_initialization.h"
//#include "elasticity/structural_element_base.h"
//#include "level_set/level_set_nonlinear_implicit_assembly.h"
//#include "elasticity/structural_nonlinear_assembly.h"
//#include "elasticity/stress_output_base.h"
//#include "elasticity/structural_near_null_vector_space.h"
//#include "base/parameter.h"
//#include "base/constant_field_function.h"
//#include "base/physics_discipline_base.h"
//#include "property_cards/solid_2d_section_element_property_card.h"
//#include "property_cards/isotropic_material_property_card.h"
//#include "boundary_condition/dirichlet_boundary_condition.h"
//#include "base/nonlinear_system.h"
//#include "level_set/level_set_system_initialization.h"
//#include "level_set/level_set_discipline.h"
//#include "level_set/level_set_transient_assembly.h"
//#include "level_set/level_set_intersection.h"
//#include "level_set/sub_cell_fe.h"
//#include "base/transient_assembly.h"
//#include "solver/first_order_newmark_transient_solver.h"
//
//// libMesh includes
//#include "libmesh/mesh_generation.h"
//#include "libmesh/exodusII_io.h"
//#include "libmesh/numeric_vector.h"
////#include "libmesh/diff_solver.h"
//#include "libmesh/nonlinear_solver.h"
//
//
//
//
//class Phi:
//public MAST::FieldFunction<Real> {
//  
//public:
//    Phi(Real l1,
//        Real l2,
//        Real rmin,
//        Real rmax):
//    MAST::FieldFunction<Real>("Phi"),
//    _l1(l1),
//    _l2(l2),
//    _rmin(rmin),
//    _rmax(rmax) {
//        
//        libmesh_assert_less(rmax, _l1*.5);
//        libmesh_assert_less(rmax, _l2*.5);
//        libmesh_assert_less(rmin,   rmax);
//    }
//    virtual ~Phi() {}
//    virtual void operator()(const libMesh::Point& p,
//                            const Real t,
//                            Real& v) const {
//        
//        libmesh_assert_less_equal(t, 1);
//        v =
//        pow(p(0)-_l1*.5, 2) +
//        pow(p(1)-_l2*.5, 2) -
//        pow((_rmin + t*(_rmax-_rmin)), 2);
//    }
//protected:
//    Real
//    _l1,
//    _l2,
//    _rmin,
//    _rmax;
//};
//
////#include "plplot/plplot.h"
//
//void plot_elem(const libMesh::Elem& e) {
//    
//    unsigned int
//    n  = e.n_nodes();
//    
//    RealVectorX
//    x  = RealVectorX::Zero(n+1),
//    y  = RealVectorX::Zero(n+1);
//    
//    for (unsigned int i=0; i<n+1; i++) {
//        x(i) = e.point(i%n)(0);
//        y(i) = e.point(i%n)(1);
//    }
//    
//    //plline(n+1, x.data(), y.data());
//}
//
//
//void plot_points(const std::vector<libMesh::Point>& pts) {
//
//    unsigned int
//    n  = pts.size();
//    
//    RealVectorX
//    x  = RealVectorX::Zero(n),
//    y  = RealVectorX::Zero(n);
//    
//    for (unsigned int i=0; i<n; i++) {
//        x(i) = pts[i](0);
//        y(i) = pts[i](1);
//    }
//    
//    //char s[] = ".";
//    //plstring(n, x.data(), y.data(), s);
//    //plpoin(n, x.data(), y.data(), -1);
//}
//
//
//class Vel: public MAST::FieldFunction<RealVectorX> {
//public:
//    Vel(): MAST::FieldFunction<RealVectorX>("vel") {}
//    
//    virtual void operator() (const libMesh::Point& p,
//                             const Real t,
//                             RealVectorX& v) const {
//        
//        v.setZero();
//        v(0) = 1.; // constant x-velocity
//    }
//    
//protected:
//    
//};
//
//
//
//MAST::PlateBendingLevelSet::PlateBendingLevelSet():
//_initialized(false),
//_length(0.),
//_width(0.),
//_mesh(nullptr),
//_eq_sys(nullptr),
//_str_sys(nullptr),
//_phi_sys(nullptr),
//_structural_sys(nullptr),
//_phi_sys_init(nullptr),
//_str_discipline(nullptr),
//_phi_discipline(nullptr),
//_th(nullptr),
//_E(nullptr),
//_nu(nullptr),
//_kappa(nullptr),
//_press(nullptr),
//_zero(nullptr),
//_th_f(nullptr),
//_E_f(nullptr),
//_nu_f(nullptr),
//_kappa_f(nullptr),
//_hoff_f(nullptr),
//_press_f(nullptr),
//_phi_vel(nullptr),
//_m_card(nullptr),
//_p_card(nullptr),
//_dirichlet_left(nullptr),
//_dirichlet_right(nullptr),
//_dirichlet_bottom(nullptr),
//_dirichlet_top(nullptr),
//_p_load(nullptr) { }
//
//
//void
//MAST::PlateBendingLevelSet::init(libMesh::ElemType e_type,
//                                 bool if_vk) {
//    
//    
//    libmesh_assert(!_initialized);
//    
//    
//    // length of domain
//    _length     = 0.3;
//    _width      = 0.3;
//    
//    
//    // create the mesh
//    _mesh       = new libMesh::ParallelMesh(this->comm());
//    
//    // initialize the mesh with one element
//    libMesh::MeshTools::Generation::build_square(*_mesh,
//                                                 20, 20,
//                                                 0, _length,
//                                                 0, _width,
//                                                 e_type);
//
//    // create the equation system
//    _eq_sys    = new  libMesh::EquationSystems(*_mesh);
//    
//    // create the libmesh system
//    _str_sys   = &(_eq_sys->add_system<MAST::NonlinearSystem>("structural"));
//    _phi_sys   = &(_eq_sys->add_system<MAST::NonlinearSystem>("level_set"));
//    
//    // FEType to initialize the system
//    libMesh::FEType fetype (libMesh::FIRST, libMesh::LAGRANGE);
//    
//    // initialize the structural system pointers
//    _structural_sys = new MAST::StructuralSystemInitialization(*_str_sys,
//                                                               _str_sys->name(),
//                                                               fetype);
//    _str_discipline = new MAST::PhysicsDisciplineBase(*_eq_sys);
//    
//    
//    // initialize the level set pointers
//    _phi_vel        = new Vel;
//    _phi_sys_init   = new MAST::LevelSetSystemInitialization(*_phi_sys,
//                                                             _phi_sys->name(),
//                                                             fetype);
//    _phi_discipline = new MAST::LevelSetDiscipline(*_eq_sys, *_phi_vel);
//
//    
////    {
////        Phi phi(_length, _width, _length*0.5*.1, _length*0.5*0.9);
////        libMesh::MeshBase::const_element_iterator
////        e_it   = _mesh->elements_begin(),
////        e_end  = _mesh->elements_end();
////
////        plsdev("xwin");
////        plinit();
////        plenv(0,.3,0,.3,0,0);
////
////        // plot the level set function
////        unsigned int npts = 100;
////        RealVectorX
////        phix = RealVectorX::Zero(npts+1),
////        phiy = RealVectorX::Zero(npts+1);
////        Real
////        pi = acos(-1.);
////
////        for (unsigned int i=0; i<npts+1; i++) {
////            phix(i) = 0.1*cos((1.*(i%npts))/(1.*npts)*2*pi)+0.15;
////            phiy(i) = 0.1*sin((1.*(i%npts))/(1.*npts)*2*pi)+0.15;
////        }
////        plline(npts+1, phix.data(), phiy.data());
////        plflush();
////
////        for ( ; e_it != e_end; e_it++) {
////
////
////            const libMesh::Elem *e = *e_it;
////            plcol0(1); // red color for elements
////            plot_elem(*e);
////            plflush();
////
////            {
////                MAST::FEBase fe(*_structural_sys);
////                fe.init(*e);
////                const std::vector<Real>& JxW = fe.get_JxW();
////                std::cout << "==== original JxW: " << std::accumulate(JxW.begin(),
////                                                                 JxW.end(), 0.) << std::endl;
////            }
////
////            MAST::LevelSetIntersection intersect;
////            intersect.init(phi, *e, 0.);
////
////            // now get elements on either side and plot
////            const std::vector<const libMesh::Elem *> &
////                    elems_low = intersect.get_sub_elems_negative_phi(),
////                    elems_hi = intersect.get_sub_elems_positive_phi();
////
////            plcol0(15); // white color for sub elements
////
////            for (unsigned int i = 0; i < elems_low.size(); i++) {
////                plot_elem(*elems_low[i]);
////                plflush();
////                // create FE
////                MAST::SubCellFE fe(*_structural_sys, intersect);
////                fe.init(*elems_low[i]);
////                const std::vector<Real>& JxW = fe.get_JxW();
////                const std::vector<libMesh::Point>& xyz = fe.get_xyz();
////                std::cout << "low: JxW: " << std::accumulate(JxW.begin(),
////                                                             JxW.end(), 0.) << std::endl;
////                plot_points(xyz);
////                plflush();
////
////                elems_low[i]->print_info();
////                for (unsigned int j=0; j<xyz.size(); j++)
////                    xyz[j].print();
////           }
////
////            for (unsigned int i=0; i<elems_hi.size(); i++) {
////                plot_elem(*elems_hi[i]);
////                plflush();
////
////                // create FE
////                MAST::SubCellFE fe(*_structural_sys, intersect);
////                fe.init(*elems_hi[i]);
////                const std::vector<Real>& JxW = fe.get_JxW();
////                const std::vector<libMesh::Point>& xyz = fe.get_xyz();
////                std::cout << "hi: JxW: " << std::accumulate(JxW.begin(),
////                                                             JxW.end(), 0.) << std::endl;
////                plot_points(fe.get_xyz());
////                plflush();
////
////                elems_hi[i]->print_info();
////                for (unsigned int j=0; j<xyz.size(); j++)
////                    xyz[j].print();
////            }
////
////        }
////        plend();
////    }
//    
//    // create and add the boundary condition and loads
//    _dirichlet_bottom = new MAST::DirichletBoundaryCondition;
//    _dirichlet_right  = new MAST::DirichletBoundaryCondition;
//    _dirichlet_top    = new MAST::DirichletBoundaryCondition;
//    _dirichlet_left   = new MAST::DirichletBoundaryCondition;
//    
//    _dirichlet_bottom->init (0, _structural_sys->vars());
//    _dirichlet_right->init  (1, _structural_sys->vars());
//    _dirichlet_top->init    (2, _structural_sys->vars());
//    _dirichlet_left->init   (3, _structural_sys->vars());
//    
//    _str_discipline->add_dirichlet_bc(0, *_dirichlet_bottom);
//    _str_discipline->add_dirichlet_bc(1,  *_dirichlet_right);
//    _str_discipline->add_dirichlet_bc(2,    *_dirichlet_top);
//    _str_discipline->add_dirichlet_bc(3,   *_dirichlet_left);
//    _str_discipline->init_system_dirichlet_bc(*_str_sys);
//    
//    // initialize the equation system
//    _eq_sys->init();
//    
//    // create the property functions and add them to the
//    
//    _th              = new MAST::Parameter("th",  0.01);
//    _E               = new MAST::Parameter("E",  72.e9);
//    _nu              = new MAST::Parameter("nu",   0.3);
//    _kappa           = new MAST::Parameter("kappa",  5./6.);
//    _zero            = new MAST::Parameter("zero",  0.);
//    _press           = new MAST::Parameter( "p",  3.e7);
//    
//    
//    
//    // prepare the vector of parameters with respect to which the sensitivity
//    // needs to be benchmarked
//    _params_for_sensitivity.push_back(_E);
//    _params_for_sensitivity.push_back(_nu);
//    _params_for_sensitivity.push_back(_th);
//    
//    
//    
//    _th_f            = new MAST::ConstantFieldFunction("h",           *_th);
//    _E_f             = new MAST::ConstantFieldFunction("E",            *_E);
//    _nu_f            = new MAST::ConstantFieldFunction("nu",          *_nu);
//    _kappa_f         = new MAST::ConstantFieldFunction("kappa",    *_kappa);
//    _hoff_f          = new MAST::ConstantFieldFunction("off",       *_zero);
//    _press_f         = new MAST::ConstantFieldFunction("pressure", *_press);
//    
//    // initialize the load
//    _p_load          = new MAST::BoundaryConditionBase(MAST::SURFACE_PRESSURE);
//    _p_load->add(*_press_f);
//    _str_discipline->add_volume_load(0, *_p_load);
//    
//    // create the material property card
//    _m_card         = new MAST::IsotropicMaterialPropertyCard;
//    
//    // add the material properties to the card
//    _m_card->add(*_E_f);
//    _m_card->add(*_nu_f);
//    _m_card->add(*_kappa_f);
//    
//    // create the element property card
//    _p_card         = new MAST::Solid2DSectionElementPropertyCard;
//    
//    // add the section properties to the card
//    _p_card->add(*_th_f);
//    _p_card->add(*_hoff_f);
//    
//    // tell the section property about the material property
//    _p_card->set_material(*_m_card);
//    _p_card->set_bending_model(MAST::MINDLIN);
//    if (if_vk) _p_card->set_strain(MAST::VON_KARMAN_STRAIN);
//    
//    _str_discipline->set_property_for_subdomain(0, *_p_card);
//        
//    _initialized = true;
//}
//
//
//
//
//
//
//
//MAST::PlateBendingLevelSet::~PlateBendingLevelSet() {
//    
//    if (_initialized) {
//        
//        delete _m_card;
//        delete _p_card;
//        
//        delete _p_load;
//        delete _dirichlet_bottom;
//        delete _dirichlet_right;
//        delete _dirichlet_top;
//        delete _dirichlet_left;
//        
//        delete _th_f;
//        delete _E_f;
//        delete _nu_f;
//        delete _kappa_f;
//        delete _hoff_f;
//        delete _press_f;
//        
//        delete _th;
//        delete _E;
//        delete _nu;
//        delete _kappa;
//        delete _zero;
//        delete _press;
//        
//        
//        
//        delete _eq_sys;
//        delete _mesh;
//        
//        delete _str_discipline;
//        delete _structural_sys;
//        delete _phi_discipline;
//        delete _phi_sys_init;
//        delete _phi_vel;
//    }
//}
//
//
//
//void
//_project_phi(MAST::FieldFunction<Real>& phi,
//             libMesh::ExplicitSystem& sys) {
//    
//    // now create a function and use it for initialization
//    class SolutionFunction:
//    public libMesh::FunctionBase<Real> {
//    public:
//        SolutionFunction(MAST::FieldFunction<Real>& phi):
//        libMesh::FunctionBase<Real>(),
//        _phi(phi) { }
//        
//        virtual std::unique_ptr<libMesh::FunctionBase<Real> > clone () const {
//            libMesh::FunctionBase<Real> *rval = new SolutionFunction(_phi);
//            return std::unique_ptr<libMesh::FunctionBase<Real> >(rval);
//        }
//        
//        // this should not get called
//        virtual Real operator()
//        (const libMesh::Point& p, const Real time) {libmesh_assert(false);}
//        
//        virtual void
//        operator() (const libMesh::Point& p,
//                    const Real time,
//                    libMesh::DenseVector<Real>& output) {
//            Real v;
//            _phi(p, time, v);
//            output(0) = v;
//        }
//    protected:
//        MAST::FieldFunction<Real>& _phi;
//    };
//    
//    SolutionFunction sol_func(phi);
//    
//    sys.project_solution(&sol_func);
//}
//
//
//
//
//const libMesh::NumericVector<Real>&
//MAST::PlateBendingLevelSet::solve(bool if_write_output) {
//    
//    
//    libmesh_assert(_initialized);
//   
//    class PhiWrapper: public MAST::FieldFunction<RealVectorX> {
//    public:
//        PhiWrapper(Phi& phi):
//        MAST::FieldFunction<RealVectorX>("phi"),
//        _phi(phi) {}
//        
//        virtual void operator()(const libMesh::Point& p,
//                                const Real t,
//                                RealVectorX& v) const {
//            Real val;
//            _phi(p, t, val);
//            v(0) = val;
//        }
//
//        Phi& _phi;
//    };
//
//    Phi phi(_length, _width, _length*0.5*.2, _length*0.5*0.9);
//    PhiWrapper phi_wrap(phi);
//
//    _phi_sys_init->initialize_solution(phi_wrap);
//    
//    // create the nonlinear assembly object
//    MAST::TransientAssembly                                  phi_assembly;
//    MAST::LevelSetTransientAssemblyElemOperations            level_set_elem_ops;
//
//    // Transient solver for time integration
//    MAST::FirstOrderNewmarkTransientSolver  phi_solver;
//    
//    // now solve the system
//    phi_assembly.set_discipline_and_system(level_set_elem_ops,
//                                              *_phi_discipline,
//                                              phi_solver,
//                                              *_phi_sys_init);
//
//    // file to write the solution for visualization
//    libMesh::ExodusII_IO exodus_writer(*_mesh);
//    
//    // time solver parameters
//    unsigned int
//    t_step            = 0;
//    
//    phi_solver.dt            = 1.e-3;
//    phi_solver.beta          = 0.5;
//    
//    // set the previous state to be same as the current state to account for
//    // zero velocity as the initial condition
//    phi_solver.solution(1).zero();
//    phi_solver.solution(1).add(1., phi_solver.solution());
//    phi_solver.solution(1).close();
//    
//    
//    if (if_write_output)
//        libMesh::out << "Writing output to : output.exo" << std::endl;
//    
//    // loop over time steps
//    while (t_step <= 100) {
//        
//        libMesh::out
//        << "Time step: " << t_step
//        << " :  t = " << _phi_sys->time
//        << std::endl;
//        
//        // write the time-step
//        if (if_write_output) {
//            
//            exodus_writer.write_timestep("output.exo",
//                                         *_eq_sys,
//                                         t_step+1,
//                                         _phi_sys->time);
//        }
//        
//        phi_solver.solve();
//        
//        phi_solver.advance_time_step();
//        t_step++;
//    }
//    
//    phi_assembly.clear_discipline_and_system();
//
//    
//    
//    /*bool if_vk = (_p_card->strain_type() == MAST::VON_KARMAN_STRAIN);
//    
//    // set the number of load steps
//    unsigned int
//    n_steps = 10;
//    //if (if_vk) n_steps = 50;
//    
//    Real
//    p0      = (*_press)();
//    
//    Phi phi(_length, _width, _length*0.5*.2, _length*0.5*0.9);
//    
//    // create the nonlinear assembly object
//    MAST::LevelSetNonlinearImplicitAssembly           assembly;
//    MAST::StructuralNonlinearAssemblyElemOperations   elem_ops;
//    
//    assembly.set_discipline_and_system(elem_ops,
//                                          *_discipline,
//                                          *_structural_sys,
//                                          phi);
//    
//    MAST::NonlinearSystem& nonlin_sys = assembly.system();
//    
//    // zero the solution before solving
//    nonlin_sys.solution->zero();
//    this->clear_stresss();
//    
//    MAST::StructuralNearNullVectorSpace nsp;
//    nonlin_sys.nonlinear_solver->nearnullspace_object = &nsp;
//    
//    
//    libMesh::ExodusII_IO exodus_writer(*_mesh);
//    
//    // now iterate over the load steps
//    for (unsigned int i=0; i<n_steps; i++) {
//        
//        //(*_press)()  =  p0*(i+1.)/(1.*n_steps);
//        nonlin_sys.time = (i+1.)/(1.*n_steps);
//        _phi_sys->time  = (i+1.)/(1.*n_steps);
//        
//        // reinit constraints for the new level set
//        nonlin_sys.reinit_constraints();
//        _project_phi(phi, *_phi_sys);
//        
//        libMesh::out
//        << "Load step: " << i << "  : t = " << nonlin_sys.time << std::endl;
//        
//        nonlin_sys.solve();
//        
//        // evaluate the outputs
//        this->clear_stresss();
//        assembly.calculate_outputs(*(_sys->solution));
//        
//        
//        if (if_write_output) {
//            
//            libMesh::out << "Writing output to : output.exo" << std::endl;
//            
//            // write the solution for visualization
//            _discipline->update_stress_strain_data();
//            exodus_writer.write_timestep("output.exo",
//                                         *_eq_sys,
//                                         i+1,
//                                         (1.*i)/(1.*(n_steps-1)));
//        }
//    }
//    assembly.clear_discipline_and_system();
//    */
//    
//    return *(_str_sys->solution);
//}
//
//
//
//
//
//const libMesh::NumericVector<Real>&
//MAST::PlateBendingLevelSet::sensitivity_solve(MAST::Parameter& p,
//                                              bool if_write_output) {
//    
//    libmesh_assert(_initialized);
//
//    Phi phi(_length, _width, _length*0.5*.1, _length*0.5*0.9);
//
//    // create the nonlinear assembly object
//    MAST::LevelSetNonlinearImplicitAssembly           assembly;
//    MAST::StructuralNonlinearAssemblyElemOperations   elem_ops;
//    
//    assembly.set_discipline_and_system(elem_ops,
//                                          *_str_discipline,
//                                          *_structural_sys,
//                                          phi);
//
//    MAST::NonlinearSystem& nonlin_sys = assembly.system();
//    
//    libMesh::ParameterVector params;
//    params.resize(1);
//    params[0]  =  p.ptr();
//    
//    // zero the solution before solving
//    nonlin_sys.add_sensitivity_solution(0).zero();
//    this->clear_stresss();
//    
//    nonlin_sys.sensitivity_solve(params);
//    
//    // evaluate sensitivity of the outputs
//    assembly.calculate_output_sensitivity(params,
//                                          true,    // true for total sensitivity
//                                          *(_str_sys->solution));
//    
//    
//    assembly.clear_discipline_and_system();
//    _str
//    
//    // write the solution for visualization
//    if (if_write_output) {
//        
//        std::ostringstream oss1, oss2;
//        oss1 << "output_"        << p.name() << ".exo";
//        oss2 << "stress_output_" << p.name() << ".exo";
//        
//        libMesh::out
//        << "Writing sensitivity output to : " << oss1.str()
//        << "  and stress/strain sensitivity to : " << oss2.str()
//        << std::endl;
//        
//        
//        _str_sys->solution->swap(_str_sys->get_sensitivity_solution(0));
//        
//        // write the solution for visualization
//        _str_discipline->update_stress_strain_data( &p);
//        libMesh::ExodusII_IO(*_mesh).write_equation_systems(oss1.str(),
//                                                            *_eq_sys);
//        
//        _str_sys->solution->swap(_str_sys->get_sensitivity_solution(0));
//    }
//    
//    return _str_sys->get_sensitivity_solution(0);
//}
//
//
