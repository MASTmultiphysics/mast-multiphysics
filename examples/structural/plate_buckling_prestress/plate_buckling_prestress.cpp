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
#include "examples/structural/plate_buckling_prestress/plate_buckling_prestress.h"
#include "elasticity/structural_system_initialization.h"
#include "elasticity/structural_element_base.h"
#include "elasticity/structural_nonlinear_assembly.h"
#include "elasticity/structural_discipline.h"
#include "elasticity/structural_system.h"
#include "elasticity/structural_buckling_eigenproblem_assembly.h"
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



/*!
 *   class defines the prestress matrix
 */
namespace MAST {
    
    class PrestressMatrix:
    public MAST::FieldFunction<RealMatrixX> {
        
    public:
        
        PrestressMatrix(const std::string& nm,
                        const MAST::Parameter& s,
                        const MAST::Parameter& p):
        MAST::FieldFunction<RealMatrixX>(nm),
        _stress(s),
        _load_param(p) {
            
            _functions.insert(&s);
            _functions.insert(&p);
        }
        
        
        virtual ~PrestressMatrix() { }
        
        
        /*!
         *    calculates the value of the function at the specified point,
         *    \par p, and time, \par t, and returns it in \p v.
         */
        virtual void operator() (const libMesh::Point& p,
                                 const Real t,
                                 RealMatrixX& v) const {
            
            v = RealMatrixX::Zero(3, 3);
            v(0,0) = _stress() * _load_param();
        }
        
        
        /*!
         *    calculates the value of the function at the specified point,
         *    \par p, and time, \par t, and returns it in \p v.
         */
        virtual void derivative (const MAST::FunctionBase& f,
                                 const libMesh::Point& p,
                                 const Real t,
                                 RealMatrixX& v) const {
            
            v = RealMatrixX::Zero(3,3);
            Real dp = 0.;
            
            // if the sensitivity parameter is the load parameter itself,
            // then the sensitivity parameter will be nonzero.
            if (_load_param.depends_on(f)) dp = 1.;
            
            v(0,0) = _stress() * dp;
        }
        
        
        
    protected:
        
        /*!
         *   parameter which defines the stress value
         */
        const MAST::Parameter& _stress;
        
        
        /*!
         *   parameter which defines the load multiplier
         */
        const MAST::Parameter& _load_param;
        
    };
}




MAST::PlateBucklingPrestress::PlateBucklingPrestress():
_initialized(false) {
    
}



void
MAST::PlateBucklingPrestress::init(libMesh::ElemType e_type,
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
    _dirichlet_bottom = new MAST::DirichletBoundaryCondition;
    _dirichlet_right  = new MAST::DirichletBoundaryCondition;
    _dirichlet_top    = new MAST::DirichletBoundaryCondition;
    _dirichlet_left   = new MAST::DirichletBoundaryCondition;
    
    std::vector<unsigned int> constrained_vars(3);
    constrained_vars[0] = 0;  // u
    constrained_vars[1] = 1;  // v
    constrained_vars[2] = 2;  // w

    _dirichlet_bottom->init (0, constrained_vars);
    _dirichlet_right->init  (1, constrained_vars);
    _dirichlet_top->init    (2, constrained_vars);
    _dirichlet_left->init   (3, constrained_vars);
    
    _discipline->add_dirichlet_bc(0, *_dirichlet_bottom);
    _discipline->add_dirichlet_bc(1,  *_dirichlet_right);
    _discipline->add_dirichlet_bc(2,    *_dirichlet_top);
    _discipline->add_dirichlet_bc(3,   *_dirichlet_left);
    _discipline->constrain_subdomain_dofs_for_var(0, 5);
    _discipline->init_system_dirichlet_bc(*_sys);
    
    // initialize the equation system
    _eq_sys->init();

    _sys->eigen_solver->set_position_of_spectrum(libMesh::LARGEST_MAGNITUDE);
    _sys->set_exchange_A_and_B(true);
    _sys->set_n_requested_eigenvalues(5);

    // create the property functions and add them to the
    
    _th              = new MAST::Parameter("th",     0.006);
    _E               = new MAST::Parameter("E",      72.e9);
    _nu              = new MAST::Parameter("nu",      0.33);
    _kappa           = new MAST::Parameter("kappa",  5./6.);
    _zero            = new MAST::Parameter("zero",      0.);
    _stress          = new MAST::Parameter("sigma", -1.0e4); // compressive stress
    _load_param      = new MAST::Parameter("load",      0.);

    
    
    // prepare the vector of parameters with respect to which the sensitivity
    // needs to be benchmarked
    _params_for_sensitivity.push_back(_E);
    _params_for_sensitivity.push_back(_nu);
    _params_for_sensitivity.push_back(_th);
    
    
    
    _th_f            = new MAST::ConstantFieldFunction("h",           *_th);
    _E_f             = new MAST::ConstantFieldFunction("E",            *_E);
    _nu_f            = new MAST::ConstantFieldFunction("nu",          *_nu);
    _kappa_f         = new MAST::ConstantFieldFunction("kappa",    *_kappa);
    _hoff_f          = new MAST::ConstantFieldFunction("off",       *_zero);
    _sigma_f         = new MAST::PrestressMatrix      ("prestress",
                                                       *_stress,
                                                       *_load_param);

    
    // create the material property card
    _m_card         = new MAST::IsotropicMaterialPropertyCard;
    
    // add the material properties to the card
    _m_card->add(*_E_f);
    _m_card->add(*_nu_f);
    _m_card->add(*_kappa_f);
    
    // create the element property card
    _p_card         = new MAST::Solid2DSectionElementPropertyCard;
    
    // add the section properties to the card
    _p_card->add(*_th_f);
    _p_card->add(*_hoff_f);
    _p_card->add(*_sigma_f);
    
    // tell the section property about the material property
    _p_card->set_material(*_m_card);
    if (if_vk) _p_card->set_strain(MAST::VON_KARMAN_STRAIN);
    
    _discipline->set_property_for_subdomain(0, *_p_card);
    
    _initialized = true;
}







MAST::PlateBucklingPrestress::~PlateBucklingPrestress() {
    
    if (_initialized) {
        
        delete _m_card;
        delete _p_card;
        
        delete _dirichlet_bottom;
        delete _dirichlet_right;
        delete _dirichlet_top;
        delete _dirichlet_left;
        delete _stress;
        delete _load_param;
        
        delete _th_f;
        delete _E_f;
        delete _nu_f;
        delete _kappa_f;
        delete _hoff_f;
        delete _sigma_f;
        
        delete _th;
        delete _E;
        delete _nu;
        delete _kappa;
        delete _zero;
        
        
        
        delete _eq_sys;
        delete _mesh;
        
        delete _discipline;
        delete _structural_sys;
    }
}



MAST::Parameter*
MAST::PlateBucklingPrestress::get_parameter(const std::string &nm) {
    
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
MAST::PlateBucklingPrestress::solve(bool if_write_output,
                                    std::vector<Real>* eig) {
    
    libmesh_assert(_initialized);

    // create the nonlinear assembly object
    MAST::StructuralBucklingEigenproblemAssembly   assembly;
    _sys->initialize_condensed_dofs(*_discipline);
    
    _sys->solution->zero();
    assembly.set_buckling_data(true,
                               *_load_param,
                               0., 1.,
                               *_sys->solution,
                               *_sys->solution);
    
    assembly.attach_discipline_and_system(*_discipline, *_structural_sys);
    _sys->eigenproblem_solve();
    
    
    // Get the number of converged eigen pairs.
    unsigned int
    nconv = std::min(_sys->get_n_converged_eigenvalues(),
                     _sys->get_n_requested_eigenvalues());
    if (eig)
        eig->resize(nconv);
    
    libMesh::ExodusII_IO*
    writer = nullptr;
    
    if (if_write_output)
        writer = new libMesh::ExodusII_IO(*_mesh);

    for (unsigned int i=0; i<nconv; i++) {
        
        // now write the eigenvalue
        Real
        re = 0.,
        im = 0.;
        _sys->get_eigenpair(i, re, im, *_sys->solution);
        re = assembly.critical_point_estimate_from_eigenproblem(re);
        
        // convert from the eigenproblem to the location where instability
        // occurs.
        if (eig)
            (*eig)[i] = re;
        
        libMesh::out
        << std::setw(35) << std::fixed << std::setprecision(15)
        << re << std::endl;
        
        if (if_write_output) {
            
            // We write the file in the ExodusII format.
            writer->write_timestep("modes.exo",
                                   *_eq_sys,
                                   i+1, i);
        }
    }
    
    assembly.clear_discipline_and_system();
}





void
MAST::PlateBucklingPrestress::sensitivity_solve(MAST::Parameter& p,
                                                std::vector<Real>& eig) {
    
    libmesh_assert(_initialized);

    _discipline->add_parameter(p);
    
    // Get the number of converged eigen pairs.
    unsigned int
    nconv = std::min(_sys->get_n_converged_eigenvalues(),
                     _sys->get_n_requested_eigenvalues());
    eig.resize(nconv);
    
    libMesh::ParameterVector params;
    params.resize(1);
    params[0]  =  p.ptr();
    
    // create the nonlinear assembly object
    MAST::StructuralBucklingEigenproblemAssembly   assembly;
    assembly.attach_discipline_and_system(*_discipline, *_structural_sys);
    _sys->eigenproblem_sensitivity_solve(params, eig);
    assembly.clear_discipline_and_system();
    
    _discipline->remove_parameter(p);
}




