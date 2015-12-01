/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2015  Manav Bhatia
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


#ifndef __mast__field_function_base__
#define __mast__field_function_base__

// C++ includes
#include <memory>

// MAST includes
#include "base/function_base.h"

// libMesh includes
#include "libmesh/point.h"
#include "libmesh/function_base.h"


namespace MAST {
    
    /*!
     *    This creates the base class for functions that have a saptial and
     *    temporal dependence, and provide sensitivity operations with respect
     *    to the functions and parameters.
     */
    template <typename ValType>
    class FieldFunction:
    public MAST::FunctionBase {
        
    public:
        FieldFunction(const std::string& nm):
        MAST::FunctionBase(nm, true)
        { }

        
        FieldFunction(const MAST::FieldFunction<ValType>& f):
        MAST::FunctionBase(f)
        { }
        
        /*!
         *   @returns a clone of the function
         */
        virtual std::auto_ptr<MAST::FieldFunction<ValType> > clone() const= 0;
        
        /*!
         *    calculates the value of the function at the specified point,
         *    \par p, and time, \par t, and returns it in \p v.
         */
        virtual void operator() (const libMesh::Point& p,
                                 const Real t,
                                 ValType& v) const = 0;
        
        
        /*!
         *    calculates the value of the function at the specified point,
         *    \par p, and time, \par t, and returns it in \p v.
         */
        virtual void derivative (const MAST::DerivativeType d,
                                 const MAST::FunctionBase& f,
                                 const libMesh::Point& p,
                                 const Real t,
                                 ValType& v) const = 0;
        
        
        /*!
         *   @returns a smart-pointer to a libMesh::FunctionBase<Real>
         *   object that could be passed to libMesh data-structures
         */
        virtual std::auto_ptr<libMesh::FunctionBase<Real> >
        libMesh_compatible_function() {libmesh_error();}
        
    protected:
    
    };
}

#endif // __mast__field_function_base__
