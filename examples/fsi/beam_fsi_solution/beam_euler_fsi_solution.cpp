/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2016  Manav Bhatia
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
#include "fluid/small_disturbance_pressure_function.h"
#include "fluid/flight_condition.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "base/boundary_condition_base.h"
#include "boundary_condition/flexible_surface_motion.h"
#include "aeroelasticity/frequency_function.h"
#include "elasticity/structural_system_initialization.h"
#include "elasticity/structural_discipline.h"
#include "elasticity/structural_element_base.h"
#include "elasticity/structural_nonlinear_assembly.h"
#include "elasticity/fsi_generalized_aero_force_assembly.h"
#include "elasticity/structural_near_null_vector_space.h"
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
#include "libmesh/centroid_partitioner.h"
#include "libmesh/nonlinear_solver.h"


extern libMesh::LibMeshInit* __init;




MAST::BeamEulerFSIAnalysis::BeamEulerFSIAnalysis():
_structural_comm(NULL),
_structural_mesh(NULL),
_fluid_mesh(NULL),
_structural_eq_sys(NULL),
_fluid_eq_sys(NULL),
_structural_sys(NULL),
_fluid_sys(NULL),
_structural_sys_init(NULL),
_structural_discipline(NULL),
_fluid_sys_init(NULL),
_fluid_discipline(NULL),
_flight_cond(NULL),
_far_field(NULL),
_symm_wall(NULL),
_slip_wall(NULL),
_pressure(NULL),
_motion_function(NULL),
_small_dist_pressure_function(NULL),
_omega(NULL),
_velocity(NULL),
_b_ref(NULL),
_omega_f(NULL),
_velocity_f(NULL),
_b_ref_f(NULL),
_freq_function(NULL),
_length(0.),
_k_lower(0.),
_k_upper(0.),
_n_k_divs(0),
_thy(NULL),
_thz(NULL),
_rho(NULL),
_E(NULL),
_nu(NULL),
_zero(NULL),
_mach(NULL),
_rho_air(NULL),
_gamma_air(NULL),
_thy_f(NULL),
_thz_f(NULL),
_rho_f(NULL),
_E_f(NULL),
_nu_f(NULL),
_hyoff_f(NULL),
_hzoff_f(NULL),
_mach_f(NULL),
_rho_air_f(NULL),
_gamma_air_f(NULL),
_m_card(NULL),
_p_card(NULL),
_dirichlet_left(NULL),
_dirichlet_right(NULL) {
    
    //////////////////////////////////////////////////////////////////////
    //    SETUP THE FLUID DATA
    //////////////////////////////////////////////////////////////////////

    // initialize the libMesh object
    _fluid_mesh              = new libMesh::ParallelMesh(__init->comm());
    _fluid_eq_sys            = new libMesh::EquationSystems(*_fluid_mesh);
    

    // tell the mesh to use the centroid partitioner with the
    // Y-coordinate as the means for sorting. This assumes that the lowest
    // layer of elements that interface with the panel are on rank 0.
    _fluid_mesh->partitioner().reset(new libMesh::CentroidPartitioner
                                     (libMesh::CentroidPartitioner::Y));
    
    
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
    _omega             = new MAST::Parameter("omega",       0.);
    _velocity          = new MAST::Parameter("velocity",  _flight_cond->velocity_magnitude);
    _b_ref             = new MAST::Parameter("b_ref",       1.);
    
    
    // now define the constant field functions based on this
    _omega_f           = new MAST::ConstantFieldFunction("omega",       *_omega);
    _velocity_f        = new MAST::ConstantFieldFunction("velocity", *_velocity);
    _b_ref_f           = new MAST::ConstantFieldFunction("b_ref",       *_b_ref);
    
    // initialize the frequency function
    _freq_function     = new MAST::FrequencyFunction("freq",
                                                     *_omega_f,
                                                     *_velocity_f,
                                                     *_b_ref_f);
    _freq_function->if_nondimensional(true);
    
    // tell the physics about boundary conditions
    _fluid_discipline->add_side_load(    panel_bc_id, *_slip_wall);
    _fluid_discipline->add_side_load( symmetry_bc_id, *_symm_wall);
    // all boundaries except the bottom are far-field
    for (unsigned int i=1; i<=3; i++)
        _fluid_discipline->add_side_load(              i, *_far_field);
    
    _small_dist_pressure_function =
    new MAST::SmallDisturbancePressureFunction(*_fluid_sys_init, *_flight_cond);

    
    _k_upper            = infile("k_upper",  0.75);
    _k_lower            = infile("k_lower",  0.05);
    _n_k_divs           = infile("n_k_divs",   10);
    
    
    //////////////////////////////////////////////////////////////////////
    //    SETUP THE STRUCTURAL DATA
    //////////////////////////////////////////////////////////////////////

    // create a communicator for the structural model on rank 0
    _structural_comm = new libMesh::Parallel::Communicator;
    if (__init->comm().rank() == 0)
        __init->comm().split(0, 0, *_structural_comm);
    else
        __init->comm().split(MPI_UNDEFINED, MPI_UNDEFINED, *_structural_comm);
    
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
    (*_b_ref) = _length;

    if (_structural_comm->get() != MPI_COMM_NULL){
        
        // create the mesh
        _structural_mesh       = new libMesh::SerialMesh(*_structural_comm);
        
        
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
        _motion_function   = new MAST::FlexibleSurfaceMotion("small_disturbance_motion",
                                                             *_structural_sys_init);
        _slip_wall->add(*_motion_function);
        
        
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
        _pressure    =  new MAST::BoundaryConditionBase(MAST::SMALL_DISTURBANCE_MOTION);
        _pressure->add(*_small_dist_pressure_function);
        _pressure->add(*_motion_function);
        _structural_discipline->add_volume_load(0, *_pressure);
    }
    
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
    
    delete _omega;
    delete _velocity;
    delete _b_ref;
    
    delete _omega_f;
    delete _velocity_f;
    delete _b_ref_f;
    
    delete _freq_function;
    
    delete _small_dist_pressure_function;
    
    
    if (_structural_comm->get() != MPI_COMM_NULL) {
        
        delete _motion_function;
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
        
        
        // delete the basis vectors
        if (_basis.size())
            for (unsigned int i=0; i<_basis.size(); i++)
                if (_basis[i]) delete _basis[i];
    }
    
    delete _structural_comm;
}





MAST::Parameter*
MAST::BeamEulerFSIAnalysis::get_parameter(const std::string &nm) {
    
    MAST::Parameter *rval = NULL;
    
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
    
    
    
    ////////////////////////////////////////////////////////////
    // STRUCTURAL MODAL EIGENSOLUTION
    ////////////////////////////////////////////////////////////
    
    MAST::StructuralNonlinearAssembly   structural_assembly;
    if (_structural_comm->get() != MPI_COMM_NULL) {
        
        // create the nonlinear assembly object
        structural_assembly.attach_discipline_and_system(*_structural_discipline,
                                                         *_structural_sys_init);
        MAST::StructuralNearNullVectorSpace nsp;
        _structural_sys->nonlinear_solver->nearnullspace_object = &nsp;
    }

    ///////////////////////////////////////////////////////////////////
    // FSI SOLUTION
    ///////////////////////////////////////////////////////////////////
    // clear flutter solver and set the output file
    MAST::MultiphysicsNonlinearSolverBase fsi_solver("fsi", 2);
    
    //fsi_solver.set_system_assembly(0,      fluid_assembly);
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

