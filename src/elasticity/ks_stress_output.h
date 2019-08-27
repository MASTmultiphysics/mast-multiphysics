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

#ifndef __mast__ks_stress_output__
#define __mast__ks_stress_output__

// C++ includes
#include <map>
#include <vector>

// MAST includes
#include "elasticity/stress_output_base.h"


// libMesh includes
#include "libmesh/elem.h"

namespace MAST {
    
    
    /*!
     *   This implements the computation of KS-constraint aggregation
     *   functional for the stress constraint.
     */
    class KSStressStrainOutput:
    public MAST::StressStrainOutputBase {
        
    public:
        
        
        /*!
         *    default constructor
         */
        KSStressStrainOutput();
        
        virtual ~KSStressStrainOutput();
        
        
        /*!
         *   calculates and returns the von Mises p-norm functional for
         *   all the elements that this object currently stores data for.
         *   This is defined as
         *   \f[  \left( \frac{\int_\Omega (\sigma_{VM}(\Omega))^p ~
         *    d\Omega}{\int_\Omega ~ d\Omega} \right)^{\frac{1}{p}} \f]
         */
        virtual void functional_for_all_elems();
        
                
        /*!
         *   calculates and returns the sensitivity of von Mises p-norm
         *   functional for the element \p e.
         */
        virtual void functional_sensitivity_for_elem
        (const MAST::FunctionBase& f,
         const libMesh::dof_id_type e_id,
         Real& dsigma_vm_val_df) const;
        
        
        /*!
         *   calculates and returns the boundary sensitivity of von Mises p-norm
         *   functional for the element \p e.
         */
        virtual void functional_boundary_sensitivity_for_elem
        (const MAST::FunctionBase& f,
         const libMesh::dof_id_type e_id,
         Real& dsigma_vm_val_df) const;
        
        
        
        /*!
         *   calculates and returns the derivative of von Mises p-norm
         *   functional wrt state vector for the specified element. This
         *   assumes that the \p von_Mises_p_norm_functional_for_all_elems()
         *   has been called to calculate the primal data.
         *   This is defined as
         *   \f[  \frac{ \frac{1}{p} \left( \int_\Omega (\sigma_{VM}(\Omega))^p ~
         *    d\Omega \right)^{\frac{1}{p}-1}}{\left(  \int_\Omega ~ d\Omega \right)^{\frac{1}{p}}}
         *    \int_\Omega p (\sigma_{VM}(\Omega))^{p-1} \frac{d \sigma_{VM}(\Omega)}{dX} ~
         *    d\Omega \f]
         */
        virtual void functional_state_derivartive_for_elem
        (const libMesh::dof_id_type e_id,
         RealVectorX& dq_dX) const;
        
        
        
    protected:
        
    };
}

#endif // __mast__ks_stress_output__
