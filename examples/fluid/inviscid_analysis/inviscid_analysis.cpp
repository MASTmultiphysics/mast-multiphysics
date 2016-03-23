
// C++ includes
#include <iostream>



// MAST includes
#include "examples/fluid/inviscid_analysis/inviscid_analysis.h"
#include "examples/fluid/meshing/panel_mesh_2D.h"
#include "fluid/conservative_fluid_system_initialization.h"
#include "fluid/conservative_fluid_discipline.h"
#include "fluid/conservative_fluid_transient_assembly.h"
#include "solver/first_order_newmark_transient_solver.h"
#include "fluid/flight_condition.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "boundary_condition/dirichlet_boundary_condition.h"

// libMesh includes
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/parameter_vector.h"
#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/elem.h"



extern libMesh::LibMeshInit* __init;



MAST::InviscidAnalysis::InviscidAnalysis() {


    // initialize the libMesh object
    _mesh              = new libMesh::ParallelMesh(__init->comm());
    _eq_sys            = new libMesh::EquationSystems(*_mesh);
    
    // add the system to be used for analysis
    _sys = &(_eq_sys->add_system<libMesh::NonlinearImplicitSystem>("fluid"));
    
    
    // initialize the flow conditions
    GetPot infile("input.in");

    
    const unsigned int
    dim                 = 2,
    nx_divs             = infile("nx_divs",          3),
    ny_divs             = infile("ny_divs",          1),
    n_max_bumps_x       = infile("n_max_bumps_x",    1),
    panel_bc_id         = infile("panel_bc_id",     10),
    symmetry_bc_id      = infile("symmetry_bc_id",  11);
    
    const bool
    if_cos_bump         = infile("if_cos_bump",  false);

    const Real
    t_by_c              = infile("t_by_c",         0.3);
    
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
    //z_divs           (nz_divs);
    
    std::auto_ptr<MeshInitializer::CoordinateDivisions>
    x_coord_divs    (new MeshInitializer::CoordinateDivisions),
    y_coord_divs    (new MeshInitializer::CoordinateDivisions);
    
    std::vector<MeshInitializer::CoordinateDivisions*>
    divs(dim);
    
    
    // now read in the values: x-coord
    if (nx_divs > 0) {
        
        for (unsigned int i_div=0; i_div<nx_divs+1; i_div++) {
            
            x_div_loc[i_div]        = infile("x_div_loc",   0., i_div);
            x_relative_dx[i_div]    = infile( "x_rel_dx",   0., i_div);
            
            if (i_div < nx_divs) //  this is only till nx_divs
                x_divs[i_div]       = infile( "x_div_nelem", 0, i_div);
        }
        
        divs[0] = x_coord_divs.get();
        x_coord_divs->init(nx_divs, x_div_loc, x_relative_dx, x_divs);
    }
    
    
    // now read in the values: y-coord
    if ((dim > 1) && (ny_divs > 0)) {
        
        for (unsigned int i_div=0; i_div<ny_divs+1; i_div++) {
            
            y_div_loc[i_div]     = infile("y_div_loc", 0., i_div);
            y_relative_dx[i_div] = infile( "y_rel_dx", 0., i_div);
            
            if (i_div < ny_divs) //  this is only till ny_divs
                y_divs[i_div]    = infile( "y_div_nelem",  0, i_div);
        }
        
        divs[1] = y_coord_divs.get();
        y_coord_divs->init(ny_divs, y_div_loc, y_relative_dx, y_divs);
    }
    


    // initialize the mesh
    MAST::PanelMesh2D().init(t_by_c,
                             if_cos_bump,
                             n_max_bumps_x,
                             panel_bc_id,
                             symmetry_bc_id,
                             divs,
                             *_mesh,
                             elem_type);
    
    _mesh->print_info();
    
    _discipline        = new MAST::ConservativeFluidDiscipline(*_eq_sys);
    _fluid_sys         = new MAST::ConservativeFluidSystemInitialization(*_sys,
                                                                         _sys->name(),
                                                                         libMesh::FEType(fe_order, fe_type),
                                                                         dim);
    
    
    // initialize the equation system for analysis
    _eq_sys->init();
    
    // print the information
    _eq_sys->print_info();
    
    // create the oundary conditions for slip-wall and far-field
    _far_field    = new MAST::BoundaryConditionBase(    MAST::FAR_FIELD),
    _slip_wall    = new MAST::BoundaryConditionBase(    MAST::SLIP_WALL);
    _symm_wall    = new MAST::BoundaryConditionBase(MAST::SYMMETRY_WALL);
    
    
    // tell the physics about these conditions
    _discipline->add_side_load(    panel_bc_id, *_slip_wall);
    _discipline->add_side_load( symmetry_bc_id, *_symm_wall);
    // all boundaries except the bottom are far-field
    for (unsigned int i=1; i<=3; i++)
        _discipline->add_side_load(              i, *_far_field);
    
        
    // time step control
    _max_time_steps    =   infile("max_time_steps", 1000);
    _time_step_size    =   infile("initial_dt",    1.e-2);
    
    
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
    _discipline->set_flight_condition(*_flight_cond);
}






MAST::InviscidAnalysis::~InviscidAnalysis() {
    
    delete _eq_sys;
    delete _mesh;
    
    delete _discipline;
    delete _fluid_sys;
    
    delete _far_field;
    delete _slip_wall;
    delete _symm_wall;
    
    delete _flight_cond;
}



MAST::Parameter*
MAST::InviscidAnalysis::get_parameter(const std::string &nm) {
    
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
        std::cout
        << std::endl
        << "Parameter not found by name: " << nm << std::endl
        << "Valid names are: "
        << std::endl;
        for (it = _params_for_sensitivity.begin(); it != end; it++)
            std::cout << "   " << (*it)->name() << std::endl;
        std::cout << std::endl;
    }
    
    return rval;
}



const libMesh::NumericVector<Real>&
MAST::InviscidAnalysis::solve(bool if_write_output) {
    
    // initialize the solution
    RealVectorX s = RealVectorX::Zero(4);
    s(0) = _flight_cond->rho();
    s(1) = _flight_cond->rho_u1();
    s(2) = _flight_cond->rho_u2();
    //s(3) = _flight_cond->rho_u3();
    s(3) = _flight_cond->rho_e();
    _fluid_sys->initialize_solution(s);
    
    // create the nonlinear assembly object
    MAST::ConservativeFluidTransientAssembly   assembly;
    
    // Transient solver for time integration
    MAST::FirstOrderNewmarkTransientSolver  solver;
    
    // now solve the system
    assembly.attach_discipline_and_system(*_discipline,
                                          solver,
                                          *_fluid_sys);
    
    libMesh::NonlinearImplicitSystem&      nonlin_sys   =
    dynamic_cast<libMesh::NonlinearImplicitSystem&>(_fluid_sys->system());

    
    // file to write the solution for visualization
    libMesh::ExodusII_IO exodus_writer(*_mesh);
    
    // time solver parameters
    unsigned int t_step  = 0;
    Real            tval = 0.;
    solver.dt            = _time_step_size;
    solver.beta          = 1.0;
    
    
    if (if_write_output)
        std::cout << "Writing output to : output.exo" << std::endl;

    // loop over time steps
    while (t_step < _max_time_steps) {
        
        std::cout
        << "Time step: " << t_step
        << " :  t = " << tval
        << " :  xdot-L2 = " << solver.velocity().l2_norm()
        << std::endl;

        // write the time-step
        if (if_write_output) {

            exodus_writer.write_timestep("output.exo",
                                         *_eq_sys,
                                         t_step+1,
                                         nonlin_sys.time);
        }
        
        solver.solve();
        
        solver.advance_time_step();
        
        tval  += solver.dt;
        t_step++;
    }
    
    assembly.clear_discipline_and_system();

    return *(_sys->solution);
}





const libMesh::NumericVector<Real>&
MAST::InviscidAnalysis::sensitivity_solve(MAST::Parameter& p,
                                          bool if_write_output) {
    
    /*_discipline->add_parameter(p);
    
    // create the nonlinear assembly object
    MAST::StructuralNonlinearAssembly   assembly;
    
    assembly.attach_discipline_and_system(*_discipline, *_structural_sys);
    
    libMesh::NonlinearImplicitSystem&      nonlin_sys   =
    dynamic_cast<libMesh::NonlinearImplicitSystem&>(assembly.system());
    
    libMesh::ParameterVector params;
    params.resize(1);
    params[0]  =  p.ptr();
    
    // zero the solution before solving
    nonlin_sys.add_sensitivity_solution(0).zero();
    this->clear_stresss();
    
    nonlin_sys.sensitivity_solve(params);
    
    // evaluate sensitivity of the outputs
    assembly.calculate_output_sensitivity(params,
                                          true,    // true for total sensitivity
                                          *(_sys->solution));
    
    
    assembly.clear_discipline_and_system();
    _discipline->remove_parameter(p);
    
    // write the solution for visualization
    if (if_write_output) {
        
        std::ostringstream oss1, oss2;
        oss1 << "output_" << p.name() << ".exo";
        oss2 << "output_" << p.name() << ".exo";
        
        std::cout
        << "Writing sensitivity output to : " << oss1.str()
        << "  and stress/strain sensitivity to : " << oss2.str()
        << std::endl;
        
        
        _sys->solution->swap(_sys->get_sensitivity_solution(0));
        
        // write the solution for visualization
        libMesh::ExodusII_IO(*_mesh).write_equation_systems(oss1.str(),
                                                            *_eq_sys);
        _discipline->plot_stress_strain_data<libMesh::ExodusII_IO>(oss2.str(), &p);
        
        _sys->solution->swap(_sys->get_sensitivity_solution(0));
    }
     */
    return _sys->get_sensitivity_solution(0);
}

