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


// MAST includes
#include "examples/structural/plate_thermally_stressed_modal_analysis/plate_thermally_stressed_modal_analysis.h"
#include "examples/structural/base/hat_stiffened_panel_mesh.h"
#include "examples/base/multilinear_interpolation.h"
#include "base/nonlinear_system.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "elasticity/structural_system_initialization.h"
#include "elasticity/structural_discipline.h"
#include "elasticity/structural_element_base.h"
#include "elasticity/structural_modal_eigenproblem_assembly.h"
#include "elasticity/stress_output_base.h"
#include "elasticity/structural_nonlinear_assembly.h"
#include "aeroelasticity/time_domain_flutter_solver.h"
#include "aeroelasticity/flutter_root_base.h"
#include "property_cards/solid_2d_section_element_property_card.h"
#include "property_cards/isotropic_material_property_card.h"
#include "boundary_condition/dirichlet_boundary_condition.h"
#include "solver/slepc_eigen_solver.h"


// libMesh includes
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/getpot.h"


extern libMesh::LibMeshInit* __init;


MAST::PlateThermallyStressedModalAnalysis::
PlateThermallyStressedModalAnalysis():
_initialized(false) {
    
}



void
MAST::PlateThermallyStressedModalAnalysis::
init(libMesh::ElemType e_type, bool if_vk) {
    
    
    libmesh_assert(!_initialized);
    
    GetPot infile("input.in");
    
    // length of domain
    _length     = infile("lx", 0.3),
    _width      = infile("ly", 0.3);
    
    _n_steps    = infile("n_steps", 15);
    
    // create the mesh
    _mesh       = new libMesh::SerialMesh(__init->comm());
    
    bool
    if_hat_stiffened = false;
    
    unsigned int
    n_stiff = 3; // a default value at this stage.
    
    std::string
    stiff_type   = infile("stiffener_type", "none");
    
    if_hat_stiffened = (stiff_type == "hat");
    
    
    if (!if_hat_stiffened)
        libMesh::MeshTools::Generation::build_square(*_mesh,
                                                     infile("nx_divs", 32),
                                                     infile("ny_divs", 32),
                                                     0, _length,
                                                     0, _width,
                                                     e_type);
    else {
        MAST::HatStiffenedPanelMesh panel_mesh;
        panel_mesh.init(n_stiff,        // n_stiff,
                        100,       // n_x_divs
                        20,       // n_y_divs_on_stiffeners
                        10,       // n_y_divs_between_stiffeners
                        _length,
                        _width,
                        infile("skin_dip_h_by_w",0.00),     // skin_dip_amplitude_by_panel_w,
                        infile("hat_dip_h_by_w",-0.00),    // hat_dip_amplitude_by_panel_w,
                        .2,       // stiff_w_by_panel_w
                        .5,       // hat_w_by_stiff_w
                        .05,      // hat_h_by_panel_w
                        *_mesh,
                        e_type);
    }
        
    // create the equation system
    _eq_sys    = new  libMesh::EquationSystems(*_mesh);
    
    // create the libmesh system
    _sys       = &(_eq_sys->add_system<MAST::NonlinearSystem>("structural"));
    _sys->set_eigenproblem_type(libMesh::GHEP);
    
    
    // FEType to initialize the system
    libMesh::FEType fetype (libMesh::FIRST, libMesh::LAGRANGE);
    
    // initialize the system to the right set of variables
    _structural_sys = new MAST::StructuralSystemInitialization(*_sys,
                                                               _sys->name(),
                                                               fetype);
    _discipline     = new MAST::StructuralDiscipline(*_eq_sys);
    
    
    // create and add the boundary condition and loads
    _dirichlet_bottom       = new MAST::DirichletBoundaryCondition;
    _dirichlet_right        = new MAST::DirichletBoundaryCondition;
    _dirichlet_top          = new MAST::DirichletBoundaryCondition;
    _dirichlet_left         = new MAST::DirichletBoundaryCondition;
    _dirichlet_stiff_left   = new MAST::DirichletBoundaryCondition;
    _dirichlet_stiff_right  = new MAST::DirichletBoundaryCondition;

    std::vector<unsigned int>
    constrained_vars(3),
    stiffener_constraints(1);
    constrained_vars[0] = 0;      // u
    constrained_vars[1] = 1;      // v
    constrained_vars[2] = 2;      // w
    
    _dirichlet_bottom->init      (0, _structural_sys->vars());
    _dirichlet_right->init       (1, _structural_sys->vars());
    _dirichlet_top->init         (2, _structural_sys->vars());
    _dirichlet_left->init        (3, _structural_sys->vars());
//    _dirichlet_bottom->init (0, _structural_sys->vars());
//    _dirichlet_right->init  (1, _structural_sys->vars());
//    _dirichlet_top->init    (2, _structural_sys->vars());
//    _dirichlet_left->init   (3, _structural_sys->vars());
    
    _discipline->add_dirichlet_bc(0,        *_dirichlet_bottom);
    _discipline->add_dirichlet_bc(1,         *_dirichlet_right);
    _discipline->add_dirichlet_bc(2,           *_dirichlet_top);
    _discipline->add_dirichlet_bc(3,          *_dirichlet_left);

    if (if_hat_stiffened) {
        
        stiffener_constraints[0] = 0; // u
        _dirichlet_stiff_left->init  (5,   stiffener_constraints);
        _dirichlet_stiff_right->init (7,   stiffener_constraints);
        _discipline->add_dirichlet_bc(5,    *_dirichlet_stiff_left);
        _discipline->add_dirichlet_bc(7,   *_dirichlet_stiff_right);
    }
    else {
        
        _dirichlet_stiff_left  = nullptr;
        _dirichlet_stiff_right = nullptr;
    }
    _discipline->init_system_dirichlet_bc(*_sys);
    
    // initialize the equation system
    _eq_sys->init();
    
    _sys->eigen_solver->set_position_of_spectrum(libMesh::LARGEST_MAGNITUDE);
    _sys->set_exchange_A_and_B(true);
    _sys->set_n_requested_eigenvalues(infile("n_eig_req", 20));
    
    
    // create the property functions and add them to the
    
    _th              = new MAST::Parameter("th",           infile("th", 0.006));
    _E               = new MAST::Parameter("E",            infile("E",  72.e9));
    _rho             = new MAST::Parameter("rho",          infile("rho",2.8e3));
    _nu              = new MAST::Parameter("nu",           infile("nu",  0.33));
    _alpha           = new MAST::Parameter("alpha",    infile("alpha", 2.5e-5));
    _kappa           = new MAST::Parameter("kappa",     infile("kappa", 5./6.));
    _zero            = new MAST::Parameter("zero",                          0.);
    _temp            = new MAST::Parameter("temperature",     infile("T", 30.));
    _press           = new MAST::Parameter( "p",              infile("p",  0.));

    
    
    // prepare the vector of parameters with respect to which the sensitivity
    // needs to be benchmarked
    _params_for_sensitivity.push_back( _E);
    _params_for_sensitivity.push_back(_nu);
    _params_for_sensitivity.push_back(_th);
    
    
    
    _th_f            = new MAST::ConstantFieldFunction("h",           *_th);
    _E_f             = new MAST::ConstantFieldFunction("E",            *_E);
    _nu_f            = new MAST::ConstantFieldFunction("nu",          *_nu);
    _rho_f           = new MAST::ConstantFieldFunction("rho",        *_rho);
    _alpha_f         = new MAST::ConstantFieldFunction("alpha_expansion", *_alpha);
    _kappa_f         = new MAST::ConstantFieldFunction("kappa",    *_kappa);
    _temp_f          = new MAST::ConstantFieldFunction("temperature", *_temp);
    _ref_temp_f      = new MAST::ConstantFieldFunction("ref_temperature", *_zero);
    _hoff_f          = new MAST::SectionOffset("off",
                                               *_th_f,
                                               infile("th_off", 0.));
    _press_f         = new MAST::ConstantFieldFunction("pressure", *_press);

    _p_load          = new MAST::BoundaryConditionBase(MAST::SURFACE_PRESSURE);
    _p_load->add(*_press_f);
    _discipline->add_volume_load(0, *_p_load);
    
    // initialize the load
    _T_load          = new MAST::BoundaryConditionBase(MAST::TEMPERATURE);
    _T_load->add(*_temp_f);
    _T_load->add(*_ref_temp_f);
    _discipline->add_volume_load(0, *_T_load);

    // create the material property card
    _m_card         = new MAST::IsotropicMaterialPropertyCard;
    
    // add the material properties to the card
    _m_card->add(*_E_f);
    _m_card->add(*_nu_f);
    _m_card->add(*_kappa_f);
    _m_card->add(*_alpha_f);
    _m_card->add(*_rho_f);
    
    
    // create the element property card
    _p_card         = new MAST::Solid2DSectionElementPropertyCard;
    
    // add the section properties to the card
    _p_card->add(*_th_f);
    _p_card->add(*_hoff_f);
    
    // tell the section property about the material property
    _p_card->set_material(*_m_card);
    if (if_vk) _p_card->set_strain(MAST::NONLINEAR_STRAIN);

    for (unsigned int i=0; i<n_stiff+1; i++)
        _discipline->set_property_for_subdomain(i, *_p_card);

    
    
    // create the output objects, one for each element
    libMesh::MeshBase::const_element_iterator
    e_it    = _mesh->elements_begin(),
    e_end   = _mesh->elements_end();
    
    // points where stress is evaluated
    std::vector<libMesh::Point> pts;
    
    if (e_type == libMesh::QUAD4 ||
        e_type == libMesh::QUAD8 ||
        e_type == libMesh::QUAD9) {
        
        pts.push_back(libMesh::Point(-1/sqrt(3), -1/sqrt(3), 1.)); // upper skin
        pts.push_back(libMesh::Point(-1/sqrt(3), -1/sqrt(3),-1.)); // lower skin
        pts.push_back(libMesh::Point( 1/sqrt(3), -1/sqrt(3), 1.)); // upper skin
        pts.push_back(libMesh::Point( 1/sqrt(3), -1/sqrt(3),-1.)); // lower skin
        pts.push_back(libMesh::Point( 1/sqrt(3),  1/sqrt(3), 1.)); // upper skin
        pts.push_back(libMesh::Point( 1/sqrt(3),  1/sqrt(3),-1.)); // lower skin
        pts.push_back(libMesh::Point(-1/sqrt(3),  1/sqrt(3), 1.)); // upper skin
        pts.push_back(libMesh::Point(-1/sqrt(3),  1/sqrt(3),-1.)); // lower skin
    }
    else if (e_type == libMesh::TRI3 ||
             e_type == libMesh::TRI6) {
        
        pts.push_back(libMesh::Point(1./3., 1./3., 1.)); // upper skin
        pts.push_back(libMesh::Point(1./3., 1./3.,-1.)); // lower skin
        pts.push_back(libMesh::Point(2./3., 1./3., 1.)); // upper skin
        pts.push_back(libMesh::Point(2./3., 1./3.,-1.)); // lower skin
        pts.push_back(libMesh::Point(1./3., 2./3., 1.)); // upper skin
        pts.push_back(libMesh::Point(1./3., 2./3.,-1.)); // lower skin
    }
    else
        libmesh_assert(false); // should not get here
    
    for ( ; e_it != e_end; e_it++) {
        
        MAST::StressStrainOutputBase * output = new MAST::StressStrainOutputBase;
        
        // tell the object to evaluate the data for this object only
        std::set<const libMesh::Elem*> e_set;
        e_set.insert(*e_it);
        output->set_elements_in_domain(e_set);
        output->set_points_for_evaluation(pts);
        output->set_volume_loads(_discipline->volume_loads());
        _outputs.push_back(output);
        
        _discipline->add_volume_output((*e_it)->subdomain_id(), *output);
    }

    
    _initialized = true;
}







MAST::PlateThermallyStressedModalAnalysis::
~PlateThermallyStressedModalAnalysis() {
    
    if (_initialized) {
        
        delete _m_card;
        delete _p_card;
        
        delete _dirichlet_bottom;
        delete _dirichlet_right;
        delete _dirichlet_top;
        delete _dirichlet_left;
        
        if ( _dirichlet_stiff_left) delete _dirichlet_stiff_left;
        if (_dirichlet_stiff_right) delete _dirichlet_stiff_right;

        delete _th_f;
        delete _E_f;
        delete _nu_f;
        delete _rho_f;
        delete _alpha_f;
        delete _kappa_f;
        delete _hoff_f;
        delete _temp_f;
        delete _ref_temp_f;

        delete _th;
        delete _E;
        delete _nu;
        delete _rho;
        delete _alpha;
        delete _kappa;
        delete _zero;
        delete _temp;
        
        delete _p_load;
        delete _press;
        delete _press_f;
        
        delete _eq_sys;
        delete _mesh;
        
        delete _discipline;
        delete _structural_sys;
        
        // iterate over the output quantities and delete them
        std::vector<MAST::StressStrainOutputBase*>::iterator
        it   =   _outputs.begin(),
        end  =   _outputs.end();
        
        for ( ; it != end; it++)
            delete *it;
        
        _outputs.clear();
    }
}



MAST::Parameter*
MAST::PlateThermallyStressedModalAnalysis::
get_parameter(const std::string &nm) {
    
    libmesh_assert(_initialized);
    
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



void
MAST::PlateThermallyStressedModalAnalysis::solve(bool if_write_output,
                                                 std::vector<Real>* eig) {
    
    libmesh_assert(_initialized);
    
    bool if_vk = (_p_card->strain_type() == MAST::NONLINEAR_STRAIN);
    
    ///////////////////////////////////////////////////////////////////
    // solve for the static deformation
    //////////////////////////////////////////////////////////////////
    
    // set the number of load steps
    unsigned int
    n_steps   = 1,
    n_perturb = 0,
    n_eig_req = _sys->get_n_requested_eigenvalues();
    if (if_vk) n_steps = _n_steps;
    
    Real
    p0      = (*_press)(),
    T0      = (*_temp)();
    
    std::ofstream freq;
    freq.open("freq.txt", std::ofstream::out);
    
    // write the header to the freq file
    freq
    << std::setw(10) << "Step"
    << std::setw(35) << "Temp"
    << std::setw(35) << "Press";
    for (unsigned int i=0; i<n_eig_req; i++) {
        std::ostringstream nm;
        nm << "mode_" << i;
        freq << std::setw(35) << nm.str();
    }
    
    
    // create the nonlinear assembly object
    MAST::StructuralNonlinearAssembly   assembly;

    // create the nonlinear assembly object
    MAST::StructuralModalEigenproblemAssembly   modal_assembly;

    // writer for the equilibrium configuration
    libMesh::ExodusII_IO exodus_writer(*_mesh);
    
    // writer for the modes
    std::vector<libMesh::ExodusII_IO*>
    mode_writer(n_eig_req);
    
    for (unsigned int i=0; i<n_eig_req; i++)
        mode_writer[i] = new libMesh::ExodusII_IO(*_mesh);
    
    //  store the base solution
    libMesh::NumericVector<Real>& base_sol = _sys->add_vector("base_solution");

    // zero the solution before solving
    _sys->solution->zero();
    
    // now iterate over the load steps
    for (int i_step=0; i_step<n_steps; i_step++) {
        
        (*_temp)()   =  T0*(i_step+1.)/(1.*n_steps);
        (*_press)()  =  p0;
 
        libMesh::out
        << "Load step: " << i_step << "  Temp = " << (*_temp)() << "  p: " << (*_press)()  << std::endl;
        
        assembly.attach_discipline_and_system(*_discipline, *_structural_sys);
        _sys->solve();
        
        freq << std::endl
        << std::setw(10) << i_step
        << std::setw(35) << std::fixed << std::setprecision(15) << (*_temp)()
        << std::setw(35) << std::fixed << std::setprecision(15) << (*_press)();
        
        // evaluate the outputs
        this->clear_stresss();
        assembly.calculate_outputs(*(_sys->solution));
        
        assembly.clear_discipline_and_system();
        
        if (if_write_output) {
            
            libMesh::out << "Writing output to : output.exo" << std::endl;
            
            // write the solution for visualization
            _discipline->update_stress_strain_data();
            exodus_writer.write_timestep("output.exo",
                                         *_eq_sys,
                                         i_step+1,
                                         (1.*i_step)/(1.*(n_steps-1)));
            
        }
        
        
        ///////////////////////////////////////////////////////////////////
        // solve for the natural frequencies
        //////////////////////////////////////////////////////////////////
        base_sol = *_sys->solution;
        
        modal_assembly.set_base_solution(base_sol);
        _sys->initialize_condensed_dofs(*_discipline);
        
        modal_assembly.attach_discipline_and_system(*_discipline, *_structural_sys);
        _sys->eigenproblem_solve();
        modal_assembly.clear_discipline_and_system();
        
        // Get the number of converged eigen pairs.
        unsigned int
        nconv = std::min(_sys->get_n_converged_eigenvalues(),
                         _sys->get_n_requested_eigenvalues());
        if (eig)
            eig->resize(nconv);
        
        for (unsigned int i=0; i<nconv; i++) {
            
            // now write the eigenvalue
            Real
            re = 0.,
            im = 0.;
            _sys->get_eigenpair(i, re, im, base_sol);
            if (eig)
                (*eig)[i] = re;
            
            libMesh::out
            << std::setw(35) << std::fixed << std::setprecision(15)
            << re << std::endl;
            
            freq << std::setw(35) << std::fixed << std::setprecision(15) << re;
            
            // if any of the eigenvalues is negative, then use the mode of the
            // smallest eigenvalue to perturb the deformation shape
            if (re < 0. && n_perturb == 0) {
                libMesh::out
                << "** Negative Frequency: Perturbing panel in first mode. **"
                << std::endl;
                
                *_sys->solution = base_sol;
                _sys->solution->scale(1./_sys->solution->linfty_norm()*1.e-1);
                n_perturb++;
            }
            
            if (if_write_output) {
            
                std::ostringstream file_name;
                
                // We write the file in the ExodusII format.
                file_name
                << "output_mode_"
                << std::setw(3) << std::setfill('0') << std::right << i
                << ".exo";

                // We write the file in the ExodusII format.
                base_sol.swap(*_sys->solution);
                mode_writer[i]->write_timestep(file_name.str(),
                                               *_eq_sys,
                                               i_step+1,
                                               (1.*i_step)/(1.*(n_steps-1)));
                base_sol.swap(*_sys->solution);
            }
        }
        
        // copy the solution back to base_sol
        base_sol = *_sys->solution;
    }
    
    // delete the exodus writes for the modes
    for (unsigned int i=0; i<n_eig_req; i++)
        delete mode_writer[i];

}





void
MAST::PlateThermallyStressedModalAnalysis::sensitivity_solve(MAST::Parameter& p,
                                                             std::vector<Real>& eig) {
    
    libmesh_assert(_initialized);
    
    ////////////////////////////////////////////////////////////////
    //   static solution sensitivity
    ////////////////////////////////////////////////////////////////
    libMesh::NumericVector<Real>
    &base_sol      = _sys->get_vector("base_solution"),
    &base_sol_sens = _sys->add_sensitivity_solution(0);
    
    _discipline->add_parameter(p);
    
    // create the nonlinear assembly object
    MAST::StructuralNonlinearAssembly   assembly;
    
    assembly.attach_discipline_and_system(*_discipline, *_structural_sys);
    
    libMesh::ParameterVector params;
    params.resize(1);
    params[0]  =  p.ptr();
    
    // zero the solution before solving
    *_sys->solution = base_sol;
    base_sol_sens.zero();
    
    _sys->sensitivity_solve(params);
    
    // evaluate sensitivity of the outputs
    assembly.calculate_output_sensitivity(params,
                                          true,    // true for total sensitivity
                                          *(_sys->solution));
    
    
    assembly.clear_discipline_and_system();

    
    ////////////////////////////////////////////////////////////////
    //   modal eigenproblem sensitivity
    ////////////////////////////////////////////////////////////////
    // Get the number of converged eigen pairs.
    unsigned int
    nconv = std::min(_sys->get_n_converged_eigenvalues(),
                     _sys->get_n_requested_eigenvalues());
    eig.resize(nconv);
    
    // create the nonlinear assembly object
    MAST::StructuralModalEigenproblemAssembly   modal_assembly;
    modal_assembly.set_base_solution(base_sol);
    modal_assembly.set_base_solution(base_sol_sens, true);
    
    modal_assembly.attach_discipline_and_system(*_discipline, *_structural_sys);
    _sys->eigenproblem_sensitivity_solve(params, eig);
    modal_assembly.clear_discipline_and_system();
    
    _discipline->remove_parameter(p);
}




void
MAST::PlateThermallyStressedModalAnalysis::clear_stresss() {
    
    libmesh_assert(_initialized);
    
    // iterate over the output quantities and delete them
    std::vector<MAST::StressStrainOutputBase*>::iterator
    it   =   _outputs.begin(),
    end  =   _outputs.end();
    
    for ( ; it != end; it++)
        (*it)->clear(false);
}

