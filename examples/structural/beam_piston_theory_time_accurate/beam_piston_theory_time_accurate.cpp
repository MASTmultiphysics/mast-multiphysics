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


// MAST includes
#include "examples/structural/beam_piston_theory_time_accurate/beam_piston_theory_time_accurate.h"
#include "elasticity/structural_system_initialization.h"
#include "elasticity/structural_element_base.h"
#include "elasticity/structural_transient_assembly.h"
#include "elasticity/structural_discipline.h"
#include "elasticity/piston_theory_boundary_condition.h"
#include "elasticity/stress_output_base.h"
#include "solver/second_order_newmark_transient_solver.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "property_cards/solid_1d_section_element_property_card.h"
#include "property_cards/isotropic_material_property_card.h"
#include "boundary_condition/dirichlet_boundary_condition.h"


// libMesh includes
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/parameter_vector.h"


extern libMesh::LibMeshInit* __init;




MAST::BeamPistonTheoryTimeAccurateAnalysis::BeamPistonTheoryTimeAccurateAnalysis() {
    
    // length of domain
    _length     = 10.;
    
    
    // create the mesh
    _mesh       = new libMesh::SerialMesh(__init->comm());
    
    // initialize the mesh with one element
    libMesh::MeshTools::Generation::build_line(*_mesh, 20, 0, _length);
    _mesh->prepare_for_use();
    
    // create the equation system
    _eq_sys    = new  libMesh::EquationSystems(*_mesh);
    
    // create the libmesh system
    _sys       = &(_eq_sys->add_system<libMesh::NonlinearImplicitSystem>("structural"));
    
    // FEType to initialize the system
    libMesh::FEType fetype (libMesh::FIRST, libMesh::LAGRANGE);
    
    // initialize the system to the right set of variables
    _structural_sys = new MAST::StructuralSystemInitialization(*_sys,
                                                               _sys->name(),
                                                               fetype);
    _discipline     = new MAST::StructuralDiscipline(*_eq_sys);

    
    // create and add the boundary condition and loads
    _dirichlet_left = new MAST::DirichletBoundaryCondition;
    _dirichlet_right= new MAST::DirichletBoundaryCondition;
    _dirichlet_left->init (0, _structural_sys->vars());
    _dirichlet_right->init(1, _structural_sys->vars());
    _discipline->add_dirichlet_bc(0, *_dirichlet_left);
    _discipline->add_dirichlet_bc(1, *_dirichlet_right);
    _discipline->init_system_dirichlet_bc(*_sys);
    
    // initialize the equation system
    _eq_sys->init();
    
    // create the property functions and add them to the
    
    _thy             = new MAST::Parameter("thy",     0.06);
    _thz             = new MAST::Parameter("thz",     0.02);
    _E               = new MAST::Parameter("E",      72.e9);
    _nu              = new MAST::Parameter("nu",      0.33);
    _rho             = new MAST::Parameter("rho",   2700.0);
    _zero            = new MAST::Parameter("zero",      0.);
    _velocity        = new MAST::Parameter("V"   ,    800.);
    _mach            = new MAST::Parameter("mach",      3.);
    _rho_air         = new MAST::Parameter("rho" ,    1.05);
    _gamma_air       = new MAST::Parameter("gamma",    1.4);
    
    
    // prepare the vector of parameters with respect to which the sensitivity
    // needs to be benchmarked
    _params_for_sensitivity.push_back(_E);
    _params_for_sensitivity.push_back(_nu);
    _params_for_sensitivity.push_back(_thy);
    _params_for_sensitivity.push_back(_thz);
    
    
    
    _thy_f           = new MAST::ConstantFieldFunction("hy",          *_thy);
    _thz_f           = new MAST::ConstantFieldFunction("hz",          *_thz);
    _E_f             = new MAST::ConstantFieldFunction("E",             *_E);
    _nu_f            = new MAST::ConstantFieldFunction("nu",           *_nu);
    _rho_f           = new MAST::ConstantFieldFunction("rho",         *_rho);
    _hyoff_f         = new MAST::ConstantFieldFunction("hy_off",     *_zero);
    _hzoff_f         = new MAST::ConstantFieldFunction("hz_off",     *_zero);
    _velocity_f      = new MAST::ConstantFieldFunction("V",      *_velocity);
    _mach_f          = new MAST::ConstantFieldFunction("mach",       *_mach);
    _rho_air_f       = new MAST::ConstantFieldFunction("rho",     *_rho_air);
    _gamma_air_f     = new MAST::ConstantFieldFunction("gamma", *_gamma_air);

    
    // now initialize the piston theory boundary conditions
    RealVectorX  vel = RealVectorX::Zero(3);
    vel(0)           = 1.;  // flow along the x-axis
    _piston_bc       = new MAST::PistonTheoryBoundaryCondition(1,     // order
                                                               vel);  // vel vector
    _piston_bc->add(*_velocity_f);
    _piston_bc->add(*_mach_f);
    _piston_bc->add(*_rho_air_f);
    _piston_bc->add(*_gamma_air_f);
    _discipline->add_volume_load(0, *_piston_bc);
    
    // create the material property card
    _m_card         = new MAST::IsotropicMaterialPropertyCard;
    
    // add the material properties to the card
    _m_card->add(  *_E_f);
    _m_card->add( *_nu_f);
    _m_card->add(*_rho_f);
    
    // create the element property card
    _p_card         = new MAST::Solid1DSectionElementPropertyCard;
    
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
        
    _discipline->set_property_for_subdomain(0, *_p_card);
    
    
    // create the output objects, one for each element
    libMesh::MeshBase::const_element_iterator
    e_it    = _mesh->elements_begin(),
    e_end   = _mesh->elements_end();
    
    // points where stress is evaluated
    std::vector<libMesh::Point> pts;
    pts.push_back(libMesh::Point(-1/sqrt(3), 1., 0.)); // upper skin
    pts.push_back(libMesh::Point(-1/sqrt(3),-1., 0.)); // lower skin
    pts.push_back(libMesh::Point( 1/sqrt(3), 1., 0.)); // upper skin
    pts.push_back(libMesh::Point( 1/sqrt(3),-1., 0.)); // lower skin
    
    for ( ; e_it != e_end; e_it++) {
        
        MAST::StressStrainOutputBase * output = new MAST::StressStrainOutputBase;
        
        // tell the object to evaluate the data for this object only
        std::set<const libMesh::Elem*> e_set;
        e_set.insert(*e_it);
        output->set_elements_in_domain(e_set);
        output->set_points_for_evaluation(pts);
        _outputs.push_back(output);
        
        _discipline->add_volume_output((*e_it)->subdomain_id(), *output);
    }
}







MAST::BeamPistonTheoryTimeAccurateAnalysis::~BeamPistonTheoryTimeAccurateAnalysis() {
    
    delete _m_card;
    delete _p_card;
    
    delete _p_load;
    delete _dirichlet_left;
    delete _dirichlet_right;
    
    delete _thy_f;
    delete _thz_f;
    delete _E_f;
    delete _nu_f;
    delete _rho_f;
    delete _hyoff_f;
    delete _hzoff_f;
    delete _velocity_f;
    delete _mach_f;
    delete _rho_air_f;
    delete _gamma_air_f;
    
    delete _thy;
    delete _thz;
    delete _E;
    delete _nu;
    delete _rho;
    delete _zero;
    delete _velocity;
    delete _mach;
    delete _rho_air;
    delete _gamma_air;
    
    
    
    delete _eq_sys;
    delete _mesh;
    
    delete _discipline;
    delete _structural_sys;

    delete _piston_bc;

    // iterate over the output quantities and delete them
    std::vector<MAST::StressStrainOutputBase*>::iterator
    it   =   _outputs.begin(),
    end  =   _outputs.end();
    
    for ( ; it != end; it++)
        delete *it;
    
    _outputs.clear();
}



MAST::Parameter*
MAST::BeamPistonTheoryTimeAccurateAnalysis::get_parameter(const std::string &nm) {
    
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
MAST::BeamPistonTheoryTimeAccurateAnalysis::solve(bool if_write_output) {
    

    // create the nonlinear assembly object
    MAST::StructuralTransientAssembly   assembly;
    
    // time solver
    MAST::SecondOrderNewmarkTransientSolver solver;
    
    assembly.attach_discipline_and_system(*_discipline,
                                          solver,
                                          *_structural_sys);
    
    libMesh::NonlinearImplicitSystem&      nonlin_sys   =
    dynamic_cast<libMesh::NonlinearImplicitSystem&>(assembly.system());
    
    // zero the solution before solving
    *nonlin_sys.solution = 1.e-4;
    this->clear_stresss();
    

    // file to write the solution for visualization
    libMesh::ExodusII_IO exodus_writer(*_mesh);
    
    
    // time solver parameters
    Real
    tval     = 0.;
    
    unsigned int
    t_step            = 0,
    n_steps           = 1000;

    solver.dt         = 1.0e-3;
    
    
    
    if (if_write_output)
        std::cout << "Writing output to : output.exo" << std::endl;
    
    // loop over time steps
    while (t_step < n_steps) {
        
        std::cout
        << "Time step: " << t_step
        << " :  t = " << tval
        << " :  xdot-L2 = " << solver.velocity().l2_norm()
        << std::endl;
        
        // write the time-step
        if (if_write_output) {
            
            // evaluate the outputs
            //assembly.calculate_outputs(*(_sys->solution));
            
            exodus_writer.write_timestep("output.exo",
                                         *_eq_sys,
                                         t_step+1,
                                         nonlin_sys.time);
            
            //_discipline->plot_stress_strain_data<libMesh::ExodusII_IO>("stress_output.exo");
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
MAST::BeamPistonTheoryTimeAccurateAnalysis::sensitivity_solve(MAST::Parameter& p,
                                     bool if_write_output) {
    
    _discipline->add_parameter(p);
    
    // create the nonlinear assembly object
    MAST::StructuralTransientAssembly   assembly;
    
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
    
    return _sys->get_sensitivity_solution(0);
}



void
MAST::BeamPistonTheoryTimeAccurateAnalysis::clear_stresss() {
    
    // iterate over the output quantities and delete them
    std::vector<MAST::StressStrainOutputBase*>::iterator
    it   =   _outputs.begin(),
    end  =   _outputs.end();
    
    for ( ; it != end; it++)
        (*it)->clear(false);
}



