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
#include "examples/thermal/bar_transient/bar_transient.h"
#include "heat_conduction/heat_conduction_system_initialization.h"
#include "heat_conduction/heat_conduction_elem_base.h"
#include "heat_conduction/heat_conduction_transient_assembly.h"
#include "heat_conduction/heat_conduction_discipline.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "property_cards/solid_1d_section_element_property_card.h"
#include "property_cards/isotropic_material_property_card.h"
#include "boundary_condition/dirichlet_boundary_condition.h"
#include "solver/first_order_newmark_transient_solver.h"
#include "base/nonlinear_system.h"


// libMesh includes
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/parameter_vector.h"


extern libMesh::LibMeshInit* __init;


MAST::BarTransient::BarTransient():
_initialized(false)  { }


void
MAST::BarTransient::init(libMesh::ElemType etype, bool if_nonlin) {
    
    libmesh_assert(!_initialized);
    
    // create the mesh
    _mesh       = new libMesh::SerialMesh(__init->comm());
    
    // initialize the mesh with one element
    libMesh::MeshTools::Generation::build_line(*_mesh, 10, 0, 10);
    
    // create the equation system
    _eq_sys    = new  libMesh::EquationSystems(*_mesh);
    
    // create the libmesh system
    _sys       = &(_eq_sys->add_system<MAST::NonlinearSystem>("conduction"));
    
    // FEType to initialize the system
    libMesh::FEType fetype (libMesh::FIRST, libMesh::LAGRANGE);
    
    // initialize the system to the right set of variables
    _thermal_sys = new MAST::HeatConductionSystemInitialization(*_sys,
                                                                _sys->name(),
                                                                fetype);
    _discipline     = new MAST::HeatConductionDiscipline(*_eq_sys);

    
    // create and add the boundary condition and loads
    _dirichlet_left      = new MAST::DirichletBoundaryCondition;
    _dirichlet_right     = new MAST::DirichletBoundaryCondition;
    _dirichlet_left->init(0, _thermal_sys->vars());
    _dirichlet_right->init(1, _thermal_sys->vars());
    _discipline->add_dirichlet_bc(0, *_dirichlet_left);
    _discipline->add_dirichlet_bc(1, *_dirichlet_right);
    _discipline->init_system_dirichlet_bc(*_sys);
    
    // initialize the equation system
    _eq_sys->init();
    
    // create the property functions and add them to the
    
    _thy             = new MAST::Parameter("thy",    0.06);
    _thz             = new MAST::Parameter("thz",    0.02);
    _k               = new MAST::Parameter("k",     190.0);
    _cp              = new MAST::Parameter("cp",    864.0);
    _rho             = new MAST::Parameter("rho",  2800.0);
    _zero            = new MAST::Parameter("zero",     0.);
    _q_source        = new MAST::Parameter( "q",     20.0);
    
    
    
    // prepare the vector of parameters with respect to which the sensitivity
    // needs to be benchmarked
    _params_for_sensitivity.push_back(_k);
    _params_for_sensitivity.push_back(_cp);
    _params_for_sensitivity.push_back(_rho);
    _params_for_sensitivity.push_back(_thy);
    _params_for_sensitivity.push_back(_thz);
    
    
    
    _thy_f           = new MAST::ConstantFieldFunction("hy",     *_thy);
    _thz_f           = new MAST::ConstantFieldFunction("hz",     *_thz);
    _k_f             = new MAST::ConstantFieldFunction("k_th",     *_k);
    _cp_f            = new MAST::ConstantFieldFunction("cp",      *_cp);
    _rho_f           = new MAST::ConstantFieldFunction("rho",    *_rho);
    _hyoff_f         = new MAST::ConstantFieldFunction("hy_off", *_zero);
    _hzoff_f         = new MAST::ConstantFieldFunction("hz_off", *_zero);
    _q_source_f      = new MAST::ConstantFieldFunction("heat_source", *_q_source);
    
    // initialize the load
    _q_source_load          = new MAST::BoundaryConditionBase(MAST::HEAT_SOURCE);
    _q_source_load->add(*_q_source_f);
    _discipline->add_volume_load(0, *_q_source_load);
    
    // create the material property card
    _m_card         = new MAST::IsotropicMaterialPropertyCard;
    
    // add the material properties to the card
    _m_card->add(*_k_f);
    _m_card->add(*_cp_f);
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
    
    _initialized = true;
}







MAST::BarTransient::~BarTransient() {

    if (!_initialized)
        return;
    
    delete _m_card;
    delete _p_card;
    
    delete _q_source_load;
    delete _dirichlet_left;
    delete _dirichlet_right;
    
    delete _thy_f;
    delete _thz_f;
    delete _k_f;
    delete _cp_f;
    delete _rho_f;
    delete _hyoff_f;
    delete _hzoff_f;
    delete _q_source_f;
    
    delete _thy;
    delete _thz;
    delete _k;
    delete _cp;
    delete _rho;
    delete _zero;
    delete _q_source;
    
    
    
    
    delete _eq_sys;
    delete _mesh;
    
    delete _discipline;
    delete _thermal_sys;
}



MAST::Parameter*
MAST::BarTransient::get_parameter(const std::string &nm) {
    
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



const libMesh::NumericVector<Real>&
MAST::BarTransient::solve(bool if_write_output) {
    
    libmesh_assert(_initialized);

    // create the nonlinear assembly object
    MAST::HeatConductionTransientAssembly   assembly;
    
    // Transient solver for time integration
    MAST::FirstOrderNewmarkTransientSolver  solver;

    assembly.attach_discipline_and_system(*_discipline,
                                          solver,
                                          *_thermal_sys);
    
    MAST::NonlinearSystem& nonlin_sys = assembly.system();
    
    // zero the solution before solving
    nonlin_sys.solution->zero();

    // file to write the solution for visualization
    libMesh::ExodusII_IO exodus_writer(*_mesh);
    
    // time solver parameters
    unsigned int t_step  = 0;
    Real            tval = 0.;
    solver.dt            = 1.0e7;
    solver.beta          = 1.0;
    
    
    if (if_write_output)
        libMesh::out << "Writing output to : output.exo" << std::endl;
    
    // ask the solver to update the initial condition for d(Temp)/dt
    solver.solve_highest_derivative_and_advance_time_step();
    
    // loop over time steps
    while (t_step < 5) {
        
        libMesh::out
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
MAST::BarTransient::sensitivity_solve(MAST::Parameter& p,
                                      bool if_write_output) {
    
    libmesh_assert(_initialized);
    
    _discipline->add_parameter(p);
    
    // create the nonlinear assembly object
    MAST::HeatConductionTransientAssembly   assembly;
    
    assembly.attach_discipline_and_system(*_discipline, *_thermal_sys);

    MAST::NonlinearSystem& nonlin_sys = assembly.system();

    libMesh::ParameterVector params;
    params.resize(1);
    params[0]  =  p.ptr();

    // zero the solution before solving
    nonlin_sys.add_sensitivity_solution(0).zero();
    
    nonlin_sys.sensitivity_solve(params);
    
    assembly.clear_discipline_and_system();
    _discipline->remove_parameter(p);
    
    // write the solution for visualization
    if (if_write_output) {
        
        std::ostringstream oss1, oss2;
        oss1 << "output_" << p.name() << ".exo";
        
        libMesh::out
        << "Writing sensitivity output to : " << oss1.str()
        << std::endl;
        
        
        _sys->solution->swap(_sys->get_sensitivity_solution(0));
        
        // write the solution for visualization
        libMesh::ExodusII_IO(*_mesh).write_equation_systems(oss1.str(),
                                                            *_eq_sys);
        
        _sys->solution->swap(_sys->get_sensitivity_solution(0));
    }
    
    return _sys->get_sensitivity_solution(0);
}


