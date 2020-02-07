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

#ifndef __mast__function_set_base__
#define __mast__function_set_base__


// C++ includes
#include <vector>
#include <map>

// MAST includes
#include "base/function_base.h"

namespace MAST {
    
    /*!
     *   provides a methods to store property values
     */
    class FunctionSetBase {
        
    public:
        FunctionSetBase();
        
        /*!
         *   destructor deletes the function pointers
         */
        virtual ~FunctionSetBase();
        

        /*!
         *   checks if the card contains the specified property value
         */
        bool contains(const std::string& nm) const;

        
        /*!
         *    adds the function to this card and returns a reference to it.
         */
        void add(MAST::FunctionBase& f);
        
        
        /*!
         *   returns a constant reference to the specified function
         */
        template <typename ValType>
        const ValType&
        get(const std::string& nm) const {
            
            std::map<std::string, MAST::FunctionBase*>::const_iterator it =
            _properties.find(nm);
            
            // make sure that this funciton exists
            if(it == _properties.end())
                libmesh_error_msg("property not found for : " << nm );
            
            return dynamic_cast<const ValType&>(*(it->second));
        }

        
        /*!
         *   returns a writable reference to the specified function
         */
        template <typename ValType>
        ValType&
        get(const std::string& nm) {
            
            std::map<std::string, MAST::FunctionBase*>::iterator it =
            _properties.find(nm);
            
            // make sure that this funciton exists
            if(it == _properties.end())
                libmesh_error_msg("property not found for : " << nm );
            
            return dynamic_cast<ValType&>(*(it->second));
        }
        
        
        /*!
         *  returns true if the property card depends on the function \p f
         */
        virtual bool depends_on(const MAST::FunctionBase& f) const;

        
    protected:
        
        /*!
         *    map of the functions in this card
         */
        std::map<std::string, MAST::FunctionBase*> _properties;
    };
    
}


#endif // __mast__function_set_base__
