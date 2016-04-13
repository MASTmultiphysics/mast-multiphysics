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
#include <fstream>
#include <iomanip>

// MAST includes
#include "examples/structural/beam_thermally_stressed_modal_analysis/beam_thermally_stressed_modal_analysis.h"
#include "elasticity/structural_system_initialization.h"
#include "elasticity/structural_element_base.h"
#include "elasticity/structural_nonlinear_assembly.h"
#include "elasticity/structural_modal_eigenproblem_assembly.h"
#include "elasticity/structural_discipline.h"
#include "elasticity/stress_output_base.h"
#include "base/parameter.h"
#include "base/nonlinear_system.h"
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


MAST::BeamThermallyStressedModalAnalysis::BeamThermallyStressedModalAnalysis(bool if_nonlin) {
    
    // length of domain
    _length     = 100.;
    
    
    // create the mesh
    _mesh       = new libMesh::SerialMesh(__init->comm());
    
    // initialize the mesh with one element
    libMesh::MeshTools::Generation::build_line(*_mesh, 20, 0, _length);
    _mesh->prepare_for_use();
    
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
    _dirichlet_left = new MAST::DirichletBoundaryCondition;
    _dirichlet_right= new MAST::DirichletBoundaryCondition;
    _dirichlet_left->init (0, _structural_sys->vars());
    _dirichlet_right->init(1, _structural_sys->vars());
    _discipline->add_dirichlet_bc(0, *_dirichlet_left);
    _discipline->add_dirichlet_bc(1, *_dirichlet_right);
    _discipline->init_system_dirichlet_bc(*_sys);
    
    // initialize the equation system
    _eq_sys->init();
    
    _sys->eigen_solver->set_position_of_spectrum(libMesh::LARGEST_MAGNITUDE);
    _sys->set_exchange_A_and_B(true);
    _sys->set_n_requested_eigenvalues(8);

    // create the property functions and add them to the
    
    _thy             = new MAST::Parameter("thy", 1.);
    _thz             = new MAST::Parameter("thz", 1.);
    _rho             = new MAST::Parameter("rho",        2);
    _E               = new MAST::Parameter("E",  30e6);
    _alpha           = new MAST::Parameter("alpha",      17.3e-6);
    _nu              = new MAST::Parameter("nu",  0.33);
    _zero            = new MAST::Parameter("zero",  0.);
    _temp            = new MAST::Parameter("temperature",  0.);
    _press           = new MAST::Parameter( "p",  1e1);
    
    
    
    // prepare the vector of parameters with respect to which the sensitivity
    // needs to be benchmarked
    _params_for_sensitivity.push_back(_E);
    _params_for_sensitivity.push_back(_nu);
    _params_for_sensitivity.push_back(_thy);
    _params_for_sensitivity.push_back(_thz);
    
    
    
    _thy_f           = new MAST::ConstantFieldFunction("hy",     *_thy);
    _thz_f           = new MAST::ConstantFieldFunction("hz",     *_thz);
    _rho_f           = new MAST::ConstantFieldFunction("rho",        *_rho);
    _E_f             = new MAST::ConstantFieldFunction("E",      *_E);
    _nu_f            = new MAST::ConstantFieldFunction("nu",     *_nu);
    _hyoff_f         = new MAST::ConstantFieldFunction("hy_off", *_zero);
    _hzoff_f         = new MAST::ConstantFieldFunction("hz_off", *_zero);
    _alpha_f         = new MAST::ConstantFieldFunction("alpha_expansion", *_alpha);
    _press_f         = new MAST::ConstantFieldFunction("pressure", *_press);
    _temp_f          = new MAST::ConstantFieldFunction("temperature", *_temp);
    _ref_temp_f      = new MAST::ConstantFieldFunction("ref_temperature", *_zero);
    
    // initialize the load
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
    _m_card->add(*_alpha_f);
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
    if (if_nonlin) _p_card->set_strain(MAST::VON_KARMAN_STRAIN);

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







MAST::BeamThermallyStressedModalAnalysis::~BeamThermallyStressedModalAnalysis() {
    
    delete _m_card;
    delete _p_card;
    
    delete _p_load;
    delete _T_load;
    delete _dirichlet_left;
    delete _dirichlet_right;
    
    delete _thy_f;
    delete _thz_f;
    delete _E_f;
    delete _nu_f;
    delete _hyoff_f;
    delete _hzoff_f;
    delete _press_f;
    delete _alpha_f;
    delete _temp_f;
    delete _ref_temp_f;
    delete _rho_f;
    
    delete _thy;
    delete _thz;
    delete _E;
    delete _nu;
    delete _zero;
    delete _press;
    delete _alpha;
    delete _temp;
    delete _rho;
    
    
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



MAST::Parameter*
MAST::BeamThermallyStressedModalAnalysis::get_parameter(const std::string &nm) {
    
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



const libMesh::NumericVector<Real>&
MAST::BeamThermallyStressedModalAnalysis::solve(bool if_write_output,
                                                std::vector<Real>* eig) {
    
    
    bool if_vk = (_p_card->strain_type() == MAST::VON_KARMAN_STRAIN);
    
    ///////////////////////////////////////////////////////////////////
    // solve for the static deformation
    //////////////////////////////////////////////////////////////////
    
    // set the number of load steps
    unsigned int
    n_steps   = 1,
    n_perturb = 0;
    if (if_vk) n_steps = 10;
    
    Real
    p0      = (*_press)(),
    T0      = (*_temp)();
    
    std::ofstream freq;
    freq.open("freq.txt", std::ofstream::out);
    
    // create the nonlinear assembly object
    MAST::StructuralNonlinearAssembly   assembly;
    
    libMesh::ExodusII_IO exodus_writer(*_mesh);
    
    // zero the solution before solving
    _sys->solution->zero();
    
    // now iterate over the load steps
    for (unsigned int i_step=0; i_step<n_steps; i_step++) {
        
        (*_temp)()   =  T0*(i_step+1.)/(1.*n_steps);
        (*_press)()  =  p0*(i_step+1.)/(1.*n_steps);
        
        libMesh::out
        << "Load step: " << i_step << "  Temp = " << (*_temp)() << "  p: " << (*_press)()  << std::endl;
        
        assembly.attach_discipline_and_system(*_discipline, *_structural_sys);
        _sys->solve();
        //}
        freq << std::endl << i_step << "   " << std::setw(35) << std::fixed << std::setprecision(15) << (*_temp)() << "  ";
        
        assembly.clear_discipline_and_system();
        
        if (if_write_output) {
            
            libMesh::out << "Writing output to : output.exo" << std::endl;
            
            // write the solution for visualization
            /*std::set<std::string> names; names.insert(_sys->name());
             libMesh::ExodusII_IO(*_mesh).write_equation_systems("output.exo",
             *_eq_sys,
             &names);*/
            exodus_writer.write_timestep("output.exo",
                                         *_eq_sys,
                                         i_step+1,
                                         (1.*i_step)/(1.*(n_steps-1)));
        }
        
        //  store the base solution
        libMesh::NumericVector<Real>& base_sol = _sys->add_vector("base_solution");
        base_sol = *_sys->solution;
        
        
        ///////////////////////////////////////////////////////////////////
        // solve for the natural frequencies
        //////////////////////////////////////////////////////////////////
        
        // create the nonlinear assembly object
        MAST::StructuralModalEigenproblemAssembly   modal_assembly;
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
            
            std::ostringstream file_name;
            
            // We write the file in the ExodusII format.
            file_name << "out_step_"
            << i_step << "_mode_"
            << std::setw(3)
            << std::setfill('0')
            << std::right
            << i
            << ".exo";
            
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
                *_sys->solution = base_sol;
                _sys->solution->scale(1./_sys->solution->linfty_norm()*1.e-1);
                n_perturb++;
            }
            
            if (if_write_output) {
                
                /* libMesh::out
                 << "Writing mode " << i << " to : "
                 << file_name.str() << std::endl;*/
                
                
                // We write the file in the ExodusII format.
                base_sol.swap(*_sys->solution);
                std::set<std::string> names; names.insert(_sys->name());
                libMesh::ExodusII_IO(*_mesh).write_equation_systems(file_name.str(),
                                                                    *_eq_sys,
                                                                    &names);
                base_sol.swap(*_sys->solution);
            }
        }
    }
    
    return *(_sys->solution);
}





const libMesh::NumericVector<Real>&
MAST::BeamThermallyStressedModalAnalysis::sensitivity_solve(MAST::Parameter& p,
                                                            std::vector<Real>& eig) {
    
    _discipline->add_parameter(p);
    
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
    
    return _sys->get_sensitivity_solution(0);
}



void
MAST::BeamThermallyStressedModalAnalysis::clear_stresss() {
    
    // iterate over the output quantities and delete them
    std::vector<MAST::StressStrainOutputBase*>::iterator
    it   =   _outputs.begin(),
    end  =   _outputs.end();
    
    for ( ; it != end; it++)
        (*it)->clear(false);
}



