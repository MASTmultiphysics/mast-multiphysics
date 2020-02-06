/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2020  Manav Bhatia and MAST authors
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

#ifndef __mast__stress_assembly__
#define __mast__stress_assembly__

// MAST includes
#include "base/assembly_base.h"

// libMesh includes
#include "libmesh/nonlinear_implicit_system.h"


namespace MAST {
    
    // Forward declerations
    class StressStrainOutputBase;
    
    
    class StressAssembly:
    public MAST::AssemblyBase {
    public:
        
        
        /*!
         *   constructor associates this assembly object with the system
         */
        StressAssembly();
        
        
        /*!
         *   destructor resets the association of this assembly object with
         *   the system
         */
        virtual ~StressAssembly();
        
        
        
        /*!
         *   updates the stresses and strains for the specified solution
         *   vector \p X. Only the maximum values out of each element are
         *   updated. This will put the stress data in the System::solution
         *   vector related to stress/strain values.
         */
        virtual void
        update_stress_strain_data(MAST::StressStrainOutputBase&       ops,
                                  const libMesh::NumericVector<Real>& X);

        /*!
         *   updates the sensitivity of stresses and strains for the specified
         *   solution vector \p X and its sensitivity, \p dXdp, with respect
         *   to parameter \p p. Only the maximum values out of each element are
         *   updated.
         */
        virtual void
        update_stress_strain_sensitivity_data(MAST::StressStrainOutputBase&       ops,
                                              const libMesh::NumericVector<Real>& X,
                                              const libMesh::NumericVector<Real>& dXdp,
                                              const MAST::FunctionBase& p,
                                              libMesh::NumericVector<Real>& dsigmadp);


    protected:
                
    };
}


#endif //__mast__stress_assembly__

