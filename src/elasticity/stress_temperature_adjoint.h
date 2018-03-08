/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2018  Manav Bhatia
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

#ifndef __mast__stress_temperature_adjoint_h__
#define __mast__stress_temperature_adjoint_h__


// MAST includes
#include "elasticity/stress_output_base.h"


// libMesh includes
#include "libmesh/elem.h"

namespace MAST {
    
    
    /*!
     *  The stress and thermoelastic analysis are dependent on temperature.
     *  This class provides the contribution to the forcing function for
     *  calculation of adjoint vector for the thermal system, as a result
     *  of the stress and themroelstic analyses.
     */
    class StressTemperatureAdjoint:
    public MAST::StressStrainOutputBase {
        
    public:
        
        /*!
         *    default constructor
         */
        StressTemperatureAdjoint();
        
        virtual ~StressTemperatureAdjoint();

        void
        set_thermal_assembly(MAST::AssemblyBase& thermal_assembly);
        
        void
        set_structural_adjoint_solution(const libMesh::NumericVector<Real>& adj_sol);
        
        virtual void output_derivative_for_elem(RealVectorX& dq_dX);
        
        
    protected:

        MAST::AssemblyBase*                            _thermal_assembly;
        std::unique_ptr<libMesh::NumericVector<Real>>  _structural_adjoint;
    };
}

#endif // __mast__stress_temperature_adjoint_h__

