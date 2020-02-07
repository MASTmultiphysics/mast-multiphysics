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

// MAST includes
#include "base/function_set_base.h"



MAST::FunctionSetBase::FunctionSetBase()
{ }
        

MAST::FunctionSetBase::~FunctionSetBase()
{ }

        

bool
MAST::FunctionSetBase::contains(const std::string& nm) const {
    
    std::map<std::string, MAST::FunctionBase*>::const_iterator it =
    _properties.find(nm);
    
    // make sure that this funciton exists
    return (it != _properties.end());
}

        
        
void
MAST::FunctionSetBase::add(MAST::FunctionBase& f) {
    
    // make sure that this funciton does not already exist
    libmesh_assert_msg(!_properties.count(f.name()),
                       "Function already exists: " + f.name());
    
    bool success = _properties.insert(std::pair<std::string, MAST::FunctionBase*>
                                      (f.name(), &f)).second;
    libmesh_assert(success);
}

        

bool
MAST::FunctionSetBase::depends_on(const MAST::FunctionBase& f) const {
    
    // check with all the properties to see if any one of them is
    // dependent on the provided parameter, or is the parameter itself
    std::map<std::string, MAST::FunctionBase*>::const_iterator
    it = _properties.begin(), end = _properties.end();
    for ( ; it!=end; it++) {
        if (it->second->depends_on(f))
            return true;
    }
    
    // if it gets here, then there is no dependency
    return false;
}
