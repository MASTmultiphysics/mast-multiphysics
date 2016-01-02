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
#include "examples/structural/plate_bending_thermal_stress/plate_bending_thermal_stress.h"
#include "examples/base/multilinear_interpolation.h"
#include "examples/structural/plate_optimization/plate_optimization_base.h"
#include "elasticity/structural_system_initialization.h"
#include "elasticity/structural_element_base.h"
#include "elasticity/structural_nonlinear_assembly.h"
#include "elasticity/structural_discipline.h"
#include "elasticity/stress_output_base.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "property_cards/solid_2d_section_element_property_card.h"
#include "property_cards/isotropic_material_property_card.h"
#include "boundary_condition/dirichlet_boundary_condition.h"

// libMesh includes
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/parameter_vector.h"


extern libMesh::LibMeshInit* __init;


MAST::PlateBendingThermalStress::PlateBendingThermalStress():
_initialized(false) {
    
}


void
MAST::PlateBendingThermalStress::init(libMesh::ElemType e_type,
                                      bool if_vk) {

    libmesh_assert(!_initialized);

    // length of domain
    _length     = 0.50,
    _width      = 0.25;
    
    
    // create the mesh
    _mesh       = new libMesh::SerialMesh(__init->comm());
    
    // initialize the mesh with one element
    libMesh::MeshTools::Generation::build_square(*_mesh,
                                                 16, 16,
                                                 0, _length,
                                                 0, _width,
                                                 e_type);
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
    std::vector<unsigned int> constrained_vars(4);
    // not constraning ty, tz will keep it simply supported
    constrained_vars[0] = 0;  // u
    constrained_vars[1] = 1;  // v
    constrained_vars[2] = 2;  // w
    constrained_vars[3] = 5;  // tz

    _dirichlet_bottom = new MAST::DirichletBoundaryCondition;
    _dirichlet_right  = new MAST::DirichletBoundaryCondition;
    _dirichlet_top    = new MAST::DirichletBoundaryCondition;
    _dirichlet_left   = new MAST::DirichletBoundaryCondition;
    
    _dirichlet_bottom->init (0, constrained_vars);
    _dirichlet_right->init  (1, constrained_vars);
    _dirichlet_top->init    (2, constrained_vars);
    _dirichlet_left->init   (3, constrained_vars);
    
    _discipline->add_dirichlet_bc(0, *_dirichlet_bottom);
    _discipline->add_dirichlet_bc(1,  *_dirichlet_right);
    _discipline->add_dirichlet_bc(2,    *_dirichlet_top);
    _discipline->add_dirichlet_bc(3,   *_dirichlet_left);
    _discipline->init_system_dirichlet_bc(*_sys);
    
    // initialize the equation system
    _eq_sys->init();
    
    // create the property functions and add them to the
    
    _th              = new MAST::Parameter("th",          0.006);
    _E               = new MAST::Parameter("E",           72.e9);
    _alpha           = new MAST::Parameter("alpha",      2.5e-5);
    _kappa           = new MAST::Parameter("kappa",       5./6.);
    _nu              = new MAST::Parameter("nu",           0.33);
    _zero            = new MAST::Parameter("zero",           0.);
    _temp            = new MAST::Parameter( "temperature",  90.);
    
    
    
    // prepare the vector of parameters with respect to which the sensitivity
    // needs to be benchmarked
    _params_for_sensitivity.push_back(_E);
    _params_for_sensitivity.push_back(_nu);
    _params_for_sensitivity.push_back(_th);
    _params_for_sensitivity.push_back(_alpha);
    
    
    _th_f            = new MAST::ConstantFieldFunction("h",           *_th);
    _E_f             = new MAST::ConstantFieldFunction("E",            *_E);
    _alpha_f         = new MAST::ConstantFieldFunction("alpha_expansion", *_alpha);
    _kappa_f         = new MAST::ConstantFieldFunction("kappa",    *_kappa);
    _nu_f            = new MAST::ConstantFieldFunction("nu",          *_nu);
    _temp_f          = new MAST::ConstantFieldFunction("temperature", *_temp);
    _ref_temp_f      = new MAST::ConstantFieldFunction("ref_temperature", *_zero);
    _hoff_f          = new MAST::SectionOffset("off",
                                               _th_f->clone().release(),
                                               1.);
    
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
    _m_card->add(*_kappa_f);
    
    // create the element property card
    _p_card         = new MAST::Solid2DSectionElementPropertyCard;
    
    // add the section properties to the card
    _p_card->add(*_th_f);
    _p_card->add(*_hoff_f);
    
    // tell the section property about the material property
    _p_card->set_material(*_m_card);
    if (if_vk) _p_card->set_strain(MAST::VON_KARMAN_STRAIN);
    
    _discipline->set_property_for_subdomain(0, *_p_card);
    
    
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
        _outputs.push_back(output);
        
        _discipline->add_volume_output((*e_it)->subdomain_id(), *output);
    }
    
    
    _initialized = true;
}







MAST::PlateBendingThermalStress::~PlateBendingThermalStress() {
    
    if (_initialized) {
        
        delete _m_card;
        delete _p_card;
        
        delete _T_load;
        delete _dirichlet_bottom;
        delete _dirichlet_right;
        delete _dirichlet_top;
        delete _dirichlet_left;
        
        delete _th_f;
        delete _E_f;
        delete _alpha_f;
        delete _nu_f;
        delete _kappa_f;
        delete _hoff_f;
        delete _temp_f;
        delete _ref_temp_f;
        
        delete _th;
        delete _E;
        delete _alpha;
        delete _kappa;
        delete _nu;
        delete _zero;
        delete _temp;
        
        
        
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
MAST::PlateBendingThermalStress::get_parameter(const std::string &nm) {
    
    libmesh_assert(_initialized);
    
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
MAST::PlateBendingThermalStress::solve(bool if_write_output) {
    
    libmesh_assert(_initialized);

    bool if_vk = (_p_card->strain_type() == MAST::VON_KARMAN_STRAIN);
    
    // set the number of load steps
    unsigned int
    n_steps = 1;
    if (if_vk) n_steps = 25;
    
    Real
    T0      = (*_temp)();

    // create the nonlinear assembly object
    MAST::StructuralNonlinearAssembly   assembly;
    
    assembly.attach_discipline_and_system(*_discipline, *_structural_sys);
    
    libMesh::NonlinearImplicitSystem&      nonlin_sys   =
    dynamic_cast<libMesh::NonlinearImplicitSystem&>(assembly.system());
    
    // zero the solution before solving
    nonlin_sys.solution->zero();
    this->clear_stresss();
    
    // now iterate over the load steps
    for (unsigned int i=0; i<n_steps; i++) {
        std::cout
        << "Load step: " << i << std::endl;
        
        (*_temp)()  =  T0*(i+1.)/(1.*n_steps);
        nonlin_sys.solve();
    }
    
    // evaluate the outputs
    assembly.calculate_outputs(*(_sys->solution));

    assembly.clear_discipline_and_system();
    
    if (if_write_output) {
        
        std::cout << "Writing output to : output.exo" << std::endl;
        
        // write the solution for visualization
        libMesh::ExodusII_IO(*_mesh).write_equation_systems("output.exo",
                                                            *_eq_sys);
        
        _discipline->plot_stress_strain_data<libMesh::ExodusII_IO>("stress_output.exo");
    }
    
    return *(_sys->solution);
}





const libMesh::NumericVector<Real>&
MAST::PlateBendingThermalStress::sensitivity_solve(MAST::Parameter& p,
                                     bool if_write_output) {
    
    libmesh_assert(_initialized);

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
MAST::PlateBendingThermalStress::clear_stresss() {
    
    libmesh_assert(_initialized);

    // iterate over the output quantities and delete them
    std::vector<MAST::StressStrainOutputBase*>::iterator
    it   =   _outputs.begin(),
    end  =   _outputs.end();
    
    for ( ; it != end; it++)
        (*it)->clear(false);
}



