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


#ifndef __mast__small_disturbance_pressure_function_h__
#define __mast__small_disturbance_pressure_function_h__


// MAST includes
#include "boundary_condition/small_disturbance_pressure.h"


// libMesh includes
#include "libmesh/system.h"
#include "libmesh/mesh_function.h"


namespace MAST {
    
    // Forward declerations
    class FrequencyFunction;
    class SystemInitialization;
    class FlightCondition;
    
    
    class SmallDisturbancePressureFunction:
    public MAST::SmallDisturbancePressure {
        
    public:
        
        SmallDisturbancePressureFunction(MAST::SystemInitialization& sys,
                                         MAST::FlightCondition&      flt);
        
        
        virtual ~SmallDisturbancePressureFunction();
        
        
        /*!
         *   initiate the mesh function for this solution
         */
        void init(libMesh::NumericVector<Real>& steady_sol,
                  libMesh::NumericVector<Real>& small_dist_sol_real,
                  libMesh::NumericVector<Real>& small_dist_sol_imag);
        
        
        /*!
         *   provides a function for the definition of surface displacement,
         *   \par w, and rotation of the surface normal in \par dn_rot
         */
        virtual void
        freq_domain_pressure(const libMesh::Point& p,
                             const bool if_cp,
                             Real&                 press,
                             Complex&              dpress);
        
        
    protected:
        
        
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
        std::auto_ptr<libMesh::MeshFunction>
        _sol_function,
        _dsol_re_function,
        _dsol_im_function;
        
        /*!
         *   steady part of solution
         */
        std::auto_ptr<libMesh::NumericVector<Real> > _sol;

        /*!
         *   real part of small-disturbance solution
         */
        std::auto_ptr<libMesh::NumericVector<Real> > _dsol_real;

        /*!
         *   imag part of small-disturbance solution
         */
        std::auto_ptr<libMesh::NumericVector<Real> > _dsol_imag;

    };
}


#endif // __mast__small_disturbance_pressure_function_h__
