/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2019  Manav Bhatia
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

#ifndef __mast_check_sensitivity_h__
#define __mast_check_sensitivity_h__

// MAST includes
#include "tests/base/test_comparisons.h"
#include "elasticity/structural_system_initialization.h"
#include "elasticity/structural_discipline.h"
#include "elasticity/structural_element_base.h"
#include "elasticity/stress_output_base.h"
#include "property_cards/solid_1d_section_element_property_card.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "property_cards/isotropic_material_property_card.h"


// libMesh includes
#include "libmesh/numeric_vector.h"


namespace MAST {
    
    template <typename ValType>
    void check_sensitivity(ValType& v) {
        
        
        const Real
        delta    = 1.e-5,
        tol      = 1.e-3;
        
        // verify the sensitivity solution of this system
        RealVectorX
        sol,
        dsol,
        dsol_fd;
        
        const libMesh::NumericVector<Real>& sol_vec = v.solve();
        
        // make sure that each stress object has a single stored value
        for (unsigned int i=0; i<v._outputs.size(); i++) {
            BOOST_CHECK((v._outputs[i]->n_elem_in_storage() == 1));
        }
        
        const unsigned int
        n_dofs     = sol_vec.size(),
        n_elems    = v._mesh->n_elem();
        
        const Real
        p_val      = 2.;
        
        // store the stress values for sensitivity calculations
        // make sure that the current setup has specified one stress output
        // per element
        libmesh_assert(v._outputs.size() == n_elems);
        RealVectorX
        stress0       =  RealVectorX::Zero(n_elems),
        dstressdp     =  RealVectorX::Zero(n_elems),
        dstressdp_fd  =  RealVectorX::Zero(n_elems);
        
        for (unsigned int i=0; i<n_elems; i++) {
            // the call to all elements should actually include a single element only
            // the p-norm used is for p=2.
            stress0(i) = v._outputs[i]->von_Mises_p_norm_functional_for_all_elems(p_val);
        }
        
        
        sol      =   RealVectorX::Zero(n_dofs);
        dsol     =   RealVectorX::Zero(n_dofs);
        
        // copy the solution for later use
        for (unsigned int i=0; i<n_dofs; i++) {
            sol(i)  = sol_vec(i);
        }
        
        // now clear the stress data structures
        v.clear_stresss();
        
        // now iterate over all the parameters and calculate the analytical
        // sensitivity and compare with the numerical sensitivity
        
        Real
        p0    = 0.,
        dp    = 0.;
        
        /////////////////////////////////////////////////////////
        // now evaluate the direct sensitivity
        /////////////////////////////////////////////////////////
        
        for (unsigned int i=0; i<v._params_for_sensitivity.size(); i++ ) {
            
            MAST::Parameter& f = *v._params_for_sensitivity[i];
            
            dsol         =   RealVectorX::Zero(n_dofs);
            dsol_fd      =   RealVectorX::Zero(n_dofs);
            dstressdp    =   RealVectorX::Zero(n_elems);
            dstressdp_fd =   RealVectorX::Zero(n_elems);
            
            // calculate the analytical sensitivity
            // analysis is required at the baseline before sensitivity solution
            // and the solution has changed after the previous perturbed solution
            v.solve();
            const libMesh::NumericVector<Real>& dsol_vec = v.sensitivity_solve(f);
            
            // make sure that each stress object has a single stored value
            for (unsigned int i=0; i<v._outputs.size(); i++) {
                BOOST_CHECK((v._outputs[i]->n_elem_in_storage() == 1));
            }
            
            // copy the sensitivity solution
            for (unsigned int j=0; j<n_dofs; j++) {
                dsol(j)  = dsol_vec(j);
            }
            
            // copy the analytical sensitivity of stress values
            for (unsigned int j=0; j<n_elems; j++) {
                dstressdp(j)  =
                v._outputs[j]->von_Mises_p_norm_functional_sensitivity_for_all_elems
                (p_val, &f);
            }
            
            // now clear the stress data structures
            v.clear_stresss();
            
            // now calculate the finite difference sensitivity
            
            // identify the perturbation in the parameter
            p0           = f();
            (fabs(p0) > 0)?  dp=delta*p0 : dp=delta;
            f()         += dp;
            
            // solve at the perturbed parameter value
            const libMesh::NumericVector<Real>& sol_vec1 = v.solve();
            
            // make sure that each stress object has a single stored value
            for (unsigned int i=0; i<v._outputs.size(); i++) {
                BOOST_CHECK((v._outputs[i]->n_elem_in_storage() == 1));
            }
            
            // copy the perturbed solution
            for (unsigned int i=0; i<n_dofs; i++) {
                dsol_fd(i)  = sol_vec1(i);
            }
            
            // calculate the finite difference sensitivity for solution
            dsol_fd -= sol;
            dsol_fd /= dp;
            
            // copy the perturbed stress values
            for (unsigned int j=0; j<n_elems; j++) {
                dstressdp_fd(j)  =
                v._outputs[j]->von_Mises_p_norm_functional_for_all_elems(p_val);
            }
            
            // calculate the finite difference sensitivity for stress
            dstressdp_fd  -= stress0;
            dstressdp_fd  /= dp;
            
            // now clear the stress data structures
            v.clear_stresss();
            
            // reset the parameter value
            f()        = p0;
            
            // now compare the solution sensitivity
            BOOST_TEST_MESSAGE("  ** dX/dp (total) wrt : " << f.name() << " **");
            BOOST_CHECK(MAST::compare_vector( dsol_fd,  dsol, tol));
            
            // now compare the stress sensitivity
            BOOST_TEST_MESSAGE("  ** dvm-stress/dp (total) wrt : " << f.name() << " **");
            BOOST_CHECK(MAST::compare_vector( dstressdp_fd, dstressdp, tol));
        }
        
        v.clear_stresss();
        
        /////////////////////////////////////////////////////////
        // now evaluate the adjoint sensitivity
        /////////////////////////////////////////////////////////
        
        
    }
    
    
}


#endif // __mast_check_sensitivity_h__
