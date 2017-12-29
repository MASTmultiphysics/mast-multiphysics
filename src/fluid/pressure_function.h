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


#ifndef __mast__pressure_function_h__
#define __mast__pressure_function_h__


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
    
    
    class PressureFunction:
    public MAST::FieldFunction<Real> {
        
    public:
        
        PressureFunction(MAST::SystemInitialization& sys,
                         MAST::FlightCondition&      flt);
        
        
        virtual ~PressureFunction();
        
        
        /*!
         *   sets the mode of the pressure value return
         */
        void set_calculate_cp(bool if_cp) {
            
            _if_cp = if_cp;
        }


        /*!
         *   sets the mode of the pressure value return
         */
        void use_reference_pressure(Real ref_press) {
            
            // only if cp is not being used
            libmesh_assert(!_if_cp);
            
            _ref_pressure = ref_press;
        }

        
        /*!
         *   initiate the mesh function for this solution
         */
        void init(const libMesh::NumericVector<Real>& steady_sol,
                  const libMesh::NumericVector<Real>* small_dist_sol = nullptr);
        
        
        /*!
         *   provides the value of the pressure at the specified point and time
         */
        virtual void
        operator() (const libMesh::Point& p,
                    const Real t,
                    Real& press) const;

        
        /*!
         *   provides the pressure perturbation. The user must have initialized
         *   the perturbed solution using the appropriate init routine.
         */
        virtual void
        perturbation(const libMesh::Point& p,
                     const Real t,
                     Real& dpress) const;

        
    protected:

        /*!
         *    the function will return cp instead of pressure if this option is
         *    true.
         */
        bool _if_cp;

        /*!
         *    the function will return pressure differential with respect to
         *    reference pressue defined in the flight condition object.
         */
        Real _ref_pressure;
        
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
        _dsol_function;
        
        /*!
         *   steady part of solution
         */
        std::unique_ptr<libMesh::NumericVector<Real> > _sol;
        
        /*!
         *   small-disturbance solution
         */
        std::unique_ptr<libMesh::NumericVector<Real> > _dsol;
    };
}


#endif // __mast__pressure_function_h__
