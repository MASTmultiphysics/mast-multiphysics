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


#ifndef __mast__frequency_domain_pressure_function_h__
#define __mast__frequency_domain_pressure_function_h__


// MAST includes
#include "base/field_function_base.h"


// libMesh includes
#include "libmesh/system.h"
#include "libmesh/mesh_function.h"


namespace MAST {
    
    // Forward declerations
    class FrequencyFunction;
    class SystemInitialization;
    class FlightCondition;
    
    
    class FrequencyDomainPressureFunction:
    public MAST::FieldFunction<Complex> {
        
    public:
        
        FrequencyDomainPressureFunction(MAST::SystemInitialization& sys,
                                         MAST::FlightCondition&      flt);
        
        
        virtual ~FrequencyDomainPressureFunction();
        
        
        /*!
         *   sets the mode of the pressure value return
         */
        void set_calculate_cp(bool if_cp) {
            
            _if_cp = if_cp;
        }

        
        /*!
         *   initiate the mesh function for this solution
         */
        void init(const libMesh::NumericVector<Real>& steady_sol,
                  const libMesh::NumericVector<Real>& small_dist_sol_real,
                  const libMesh::NumericVector<Real>& small_dist_sol_imag);
        
        
        /*!
         *   provides the complex pressure perturbation
         */
        virtual void
        operator() (const libMesh::Point& p,
                    const Real t,
                    Complex& dp) const;
        
        
    protected:
        
        /*!
         *    the function will return cp instead of pressure if this option is
         *    true.
         */
        bool _if_cp;
        
        /*!
         *   system associated with the mesh and solution vector
         */
        MAST::SystemInitialization&         _system;
        
        
        /*!
         *    flight condition
         */
        MAST::FlightCondition&              _flt_cond;
        
        
        /*!
         *   mesh function that interpolates the solution
         */
        std::unique_ptr<libMesh::MeshFunction>
        _sol_function,
        _dsol_re_function,
        _dsol_im_function;
        
        /*!
         *   steady part of solution
         */
        std::unique_ptr<libMesh::NumericVector<Real> > _sol;

        /*!
         *   real part of small-disturbance solution
         */
        std::unique_ptr<libMesh::NumericVector<Real> > _dsol_real;

        /*!
         *   imag part of small-disturbance solution
         */
        std::unique_ptr<libMesh::NumericVector<Real> > _dsol_imag;

    };
}


#endif // __mast__frequency_domain_pressure_function_h__
