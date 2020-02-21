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
     *    to the functions and parameters. The field-function can provide the
     *    following:
     *       - function value,
     *       - small perturbations to support linearized analyses,
     *       - sensitivity with respect to a specified function.
     */
    template <typename ValType>
    class FieldFunction:
    public MAST::FunctionBase {
        
    public:
        FieldFunction(const std::string& nm):
        MAST::FunctionBase(nm, true)
        { }

        
        /*!
         *    calculates the value of the function and returns it in \p v.
         */
        virtual void operator() (ValType& v) const {
            
            libmesh_error_msg("Must be implemented in derived class; " << __PRETTY_FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__);
        }
        
        
        /*!
         *    calculates the perturbation and returns it in \p v.
         */
        virtual void perturbation (ValType& v) const {
            
            libmesh_error_msg("Must be implemented in derived class; " << __PRETTY_FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__);
        }

        
        /*!
         *    calculates the value of the function derivative and
         *    returns it in \p v.
         */
        virtual void derivative (const MAST::FunctionBase& f,
                                 ValType& v) const {
            
            libmesh_error_msg("Must be implemented in derived class; " << __PRETTY_FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__);
        }

        
        /*!
         *    calculates the value of the function at the specified point,
         *    \p p, and time, \p t, and returns it in \p v.
         */
        virtual void operator() (const libMesh::Point& p,
                                 const Real t,
                                 ValType& v) const {
            
            libmesh_error_msg("Must be implemented in derived class; " << __PRETTY_FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__);
        }
        
        
        /*!
         *    calculates the value of a perturbation in function at the 
         *    specified point, \p p, and time, \p t, and returns it
         *    in \p v.
         */
        virtual void perturbation (const libMesh::Point& p,
                                   const Real t,
                                   ValType& v) const {
            
            libmesh_error_msg("Must be implemented in derived class; " << __PRETTY_FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__);
        }

        
        /*!
         *    calculates the value of the derivative of function with respect to
         *    the function \p f at the specified point, \p p, and time,
         *    \p t, and returns it in \p v.
         */
        virtual void derivative (const MAST::FunctionBase& f,
                                 const libMesh::Point& p,
                                 const Real t,
                                 ValType& v) const {
            
            libmesh_error_msg("Must be implemented in derived class; " << __PRETTY_FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__);
        }
        
    protected:
    
    };
}

#endif // __mast__field_function_base__
