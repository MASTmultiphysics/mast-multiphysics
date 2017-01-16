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

#ifndef __mast_linearized_conservative_fluid_elem__
#define __mast_linearized_conservative_fluid_elem__


// MAST includes
#include "fluid/conservative_fluid_element_base.h"


namespace MAST {
    
    // Forward declerations
    class FrequencyFunction;
    
    
    class LinearizedConservativeFluidElem:
    public MAST::ConservativeFluidElementBase {
        
    public:
        
        LinearizedConservativeFluidElem(MAST::SystemInitialization& sys,
                                        const libMesh::Elem& elem,
                                        const MAST::FlightCondition& f);
        
        
        virtual ~LinearizedConservativeFluidElem();
        
        
        /*!
         *   internal force contribution to system residual
         */
        virtual bool
        internal_residual (bool request_jacobian,
                           RealVectorX& f,
                           RealMatrixX& jac);
        
        /*!
         *   side external force contribution to system residual
         */
        bool
        side_external_residual (bool request_jacobian,
                                RealVectorX& f,
                                RealMatrixX& jac,
                                std::multimap<libMesh::boundary_id_type, MAST::BoundaryConditionBase*>& bc);
        
                
        
    protected:
        
        
        /*!
         *    residual of the slip wall that may be oscillating.
         */
        virtual bool
        slip_wall_surface_residual(bool request_jacobian,
                                   RealVectorX& f,
                                   RealMatrixX& jac,
                                   const unsigned int s,
                                   MAST::BoundaryConditionBase& p);
        
        
    };
}



#endif // __mast_linearized_conservative_fluid_elem__
