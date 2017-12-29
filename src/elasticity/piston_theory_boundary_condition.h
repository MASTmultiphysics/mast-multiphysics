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

#ifndef __mast__piston_theory_boundary_condition__
#define __mast__piston_theory_boundary_condition__

// MAST includes
#include "base/boundary_condition_base.h"
#include "base/field_function_base.h"


namespace MAST {
    
    // Forward declerations
    class Parameter;
    
    // defines the function class the evaluates the pressure for calculation
    // of surface loads. The deflection and velocity data is provided
    // for calculation
    class PistonTheoryPressure:
    public MAST::FieldFunction<Real> {
        
    public:
        
        PistonTheoryPressure(unsigned int order,
                             const MAST::FieldFunction<Real>& V,
                             const MAST::FieldFunction<Real>& M,
                             const MAST::FieldFunction<Real>& rho,
                             const MAST::FieldFunction<Real>& gamma,
                             const MAST::FieldFunction<Real>& dwdx,
                             const MAST::FieldFunction<Real>& dwdt);
        
        virtual ~PistonTheoryPressure();
        
        
        virtual void operator() (const libMesh::Point& p,
                                 const Real t,
                                 Real& m) const;
        
        virtual void derivative (const MAST::FunctionBase& f,
                                 const libMesh::Point& p,
                                 const Real t,
                                 Real& m) const;
        
    protected:
        
        const unsigned int  _order;
        
        const MAST::FieldFunction<Real>
        &_V_inf,
        &_M_inf,
        &_rho_inf,
        &_gamma,
        &_dwdx,
        &_dwdt;
    };
    
    
    
    
    
    // defines the function class the evaluates the pressure for calculation
    // of surface loads. The deflection and velocity data is provided
    // for calculation
    class PistonTheoryPressureXDerivative:
    public MAST::FieldFunction<Real> {
        
    public:
        
        PistonTheoryPressureXDerivative(unsigned int order,
                                        const MAST::FieldFunction<Real>& V,
                                        const MAST::FieldFunction<Real>& M,
                                        const MAST::FieldFunction<Real>& rho,
                                        const MAST::FieldFunction<Real>& gamma,
                                        const MAST::FieldFunction<Real>& dwdx,
                                        const MAST::FieldFunction<Real>& dwdt);
        
        
        virtual ~PistonTheoryPressureXDerivative();
        
        
        virtual void operator() (const libMesh::Point& p,
                                 const Real t,
                                 Real& m) const;
        
        virtual void derivative (const MAST::FunctionBase& f,
                                 const libMesh::Point& p,
                                 const Real t,
                                 Real& m) const;
        
    protected:
        
        const unsigned int  _order;
        
        const MAST::FieldFunction<Real>
        &_V_inf,
        &_M_inf,
        &_rho_inf,
        &_gamma,
        &_dwdx,
        &_dwdt;
    };
    
    
    
    
    // defines the function class the evaluates the pressure for calculation
    // of surface loads. The deflection and velocity data is provided
    // for calculation
    class PistonTheoryPressureXdotDerivative:
    public MAST::FieldFunction<Real> {
        
    public:
        
        PistonTheoryPressureXdotDerivative(unsigned int order,
                                           const MAST::FieldFunction<Real>& V,
                                           const MAST::FieldFunction<Real>& M,
                                           const MAST::FieldFunction<Real>& rho,
                                           const MAST::FieldFunction<Real>& gamma,
                                           const MAST::FieldFunction<Real>& dwdx,
                                           const MAST::FieldFunction<Real>& dwdt);
        
        
        virtual ~PistonTheoryPressureXdotDerivative();
        
        
        virtual void operator() (const libMesh::Point& p,
                                 const Real t,
                                 Real& m) const;
        
        virtual void derivative (const MAST::FunctionBase& f,
                                 const libMesh::Point& p,
                                 const Real t,
                                 Real& m) const;
        
    protected:
        
        const unsigned int  _order;
        
        const MAST::FieldFunction<Real>
        &_V_inf,
        &_M_inf,
        &_rho_inf,
        &_gamma,
        &_dwdx,
        &_dwdt;
    };
    
    
    
    class PistonTheoryBoundaryCondition:
    public MAST::BoundaryConditionBase {
        
    public:
        
        /*!
         *  Constructor for the Piston Theory boundary condition
         *  object. The arguments needed for initialization are
         *  \p order: order of piston theory, \p mach: mach number
         *  \p a_inf: ambient speed of sound, \p gamma: ratio of
         *  specific heats at constant pressure and constant volume,
         *  \p rho: ambient density for calculation of dynamic pressure,
         *  \p vel_vec: velocity unit vector
         */
        PistonTheoryBoundaryCondition(unsigned int order,
                                      const RealVectorX& vel_vec);
        
        
        virtual ~PistonTheoryBoundaryCondition();
        
        
        /*!
         *  @returns the order of piston theory to be used
         */
        unsigned int order() const;
        
        /*!
         *  @returns velocity vector
         */
        const RealVectorX& vel_vec() const;
        
        
        /*!
         *   @returns a smart-pointer to the pressure function
         */
        std::unique_ptr<MAST::FieldFunction<Real> >
        get_pressure_function(const MAST::FieldFunction<Real>& dwdx,
                              const MAST::FieldFunction<Real>& dwdt) const;

        
        
        /*!
         *   @returns a smart-pointer to the pressure function
         */
        std::unique_ptr<MAST::FieldFunction<Real> >
        get_dpdx_function(const MAST::FieldFunction<Real>& dwdx,
                          const MAST::FieldFunction<Real>& dwdt) const;

        
        
        /*!
         *   @returns a smart-pointer to the pressure function
         */
        std::unique_ptr<MAST::FieldFunction<Real> >
        get_dpdxdot_function(const MAST::FieldFunction<Real>& dwdx,
                             const MAST::FieldFunction<Real>& dwdt) const;

        
    protected:
        
        
        /*!
         *   Order of the boundary condition
         */
        unsigned int _order;
        
        
        /*!
         *   Ambient flow velocity vector
         */
        RealVectorX _vel_vec;
        
    };
}


#endif /* __mast__piston_theory_boundary_condition__ */

