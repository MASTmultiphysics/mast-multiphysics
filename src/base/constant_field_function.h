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

#ifndef __mast__constant_field_function__
#define __mast__constant_field_function__

// MAST includes
#include "base/field_function_base.h"


namespace MAST {
    
    // Forward declerations
    class Parameter;
    
    class ConstantFieldFunction:
    public MAST::FieldFunction<Real> {
    
    public:
        
        ConstantFieldFunction(const std::string& nm,
                              const MAST::Parameter& p);
        
        
        virtual ~ConstantFieldFunction();

        
        /*!
         *    calculates the value of the function at the specified point,
         *    \par p, and time, \par t, and returns it in \p v.
         */
        virtual void operator() (Real& v) const;
        
        
        /*!
         *    calculates the value of the function at the specified point,
         *    \par p, and time, \par t, and returns it in \p v.
         */
        virtual void derivative (const MAST::DerivativeType d,
                                 const MAST::FunctionBase& f,
                                 Real& v) const;

        
        
        /*!
         *    calculates the value of the function at the specified point,
         *    \par p, and time, \par t, and returns it in \p v.
         */
        virtual void operator() (const libMesh::Point& p,
                                 const Real t,
                                 Real& v) const;
        
        
        /*!
         *    calculates the value of the function at the specified point,
         *    \par p, and time, \par t, and returns it in \p v.
         */
        virtual void derivative (const MAST::DerivativeType d,
                                 const MAST::FunctionBase& f,
                                 const libMesh::Point& p,
                                 const Real t,
                                 Real& v) const;

        
        
    protected:

        /*!
         *   parameter which defines this field function
         */
        const MAST::Parameter& _p;
        
    };
}


#endif // __mast__constant_field_function__
