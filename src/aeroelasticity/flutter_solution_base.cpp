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


// MAST includes
#include "aeroelasticity/flutter_solution_base.h"
#include "aeroelasticity/flutter_root_base.h"


MAST::FlutterSolutionBase::~FlutterSolutionBase() {
    
    std::vector<MAST::FlutterRootBase*>::iterator
    it = _roots.begin();
    for ( ; it != _roots.end(); it++)
        delete *it;
}




void
MAST::FlutterSolutionBase::swap_root(MAST::FlutterSolutionBase& sol,
                                     unsigned int root_num) {
    
    libmesh_assert(root_num < _roots.size());
    libmesh_assert(root_num < sol._roots.size());
    
    std::swap(_roots[root_num], sol._roots[root_num]);
}


