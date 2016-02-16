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

#ifndef __mast__function_base__
#define __mast__function_base__

// C++ includes
#include <set>


//  MAST includes
#include "base/mast_data_types.h"


namespace MAST
{
    
    /*!
     *   Enumeration describing the derivative type
     */
    enum DerivativeType {
        PARTIAL_DERIVATIVE,
        TOTAL_DERIVATIVE
    };
    
    
    class FunctionBase {
    public:
        
        /*!
         *   initializes the parameter to the given name
         */
        FunctionBase(const std::string& nm ,
                     const bool is_field_func):
        _name(nm),
        _is_field_func(is_field_func),
        _master(NULL)
        { }

        
        /*!
         *    Copy constructor
         */
        FunctionBase(const MAST::FunctionBase& f):
        _name(f._name),
        _is_field_func(f._is_field_func),
        _master(f.master())
        { }

        
        /*!
         *   virtual destructor
         */
        virtual ~FunctionBase() { }
        
        
        /*!
         *   returns the name of this function
         */
        const std::string& name() const {
            return _name;
        }
        
        
        const MAST::FunctionBase* master() const {
            if (_master)     // this function has a master
                return _master;
            else             // this is the master
                return this;
        }
                
        
        /*!
         *  returns true if the function depends on the provided value
         */
        virtual bool depends_on(const MAST::FunctionBase& f) const {
            if ((_functions.count(f.master())) || // one of the functions is the master
                (f.master() == this->master()))   // this function is the same
                return true;
            
            // check with all functions if they are dependent
            std::set<const MAST::FunctionBase*>::const_iterator
            it = _functions.begin(), end = _functions.end();
            
            for ( ; it != end; it++)
                if ((*it)->depends_on(*f.master()))
                    return true;
            
            // if it gets here, then there is no dependency
            return false;
        }
        
        
        /*!
         *  @returns true if the function is a shape parameter. False by
         *  default. This should be reimplemneted in a new function
         *  that is a shape function.
         */
        virtual bool is_shape_parameter() const {
            return false;
        }
        
    protected:
        
        /*!
         *    name of this parameter
         */
        std::string _name;
        
        /*!
         *    flag to store the nature of field function
         */
        bool _is_field_func;
        
        /*!
         *    pointer to the master function that this is a copy of. A NULL
         *    pointer implies that this is the master function
         */
        const MAST::FunctionBase* _master;
        
        /*!
         *   set of functions that \p this function depends on
         */
        std::set<const MAST::FunctionBase*> _functions;
    };
    
}

#endif // __mast__function_base__
