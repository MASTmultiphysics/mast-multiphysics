/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2014  Manav Bhatia
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
#include "base/function_base.h"



//bool
//MAST::FunctionBase::depends_on(const MAST::FunctionBase &p) const {
//    
//    // only first order sensitivities are calculated at this point
//    libmesh_assert_equal_to(p.total_order(), 1);
//    
//    const MAST::FunctionBase::ParameterMap& p_map = p.get_map();
//    MAST::SensitivityParameters::ParameterMap::const_iterator it, end;
//    it = p_map.begin(); end = p_map.end();
//    
//    const MAST::FunctionBase& f = *(it->first);
//    return this->depends_on(f);
//}




