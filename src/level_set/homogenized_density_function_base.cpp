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
#include "level_set/homogenized_density_function_base.h"
#include "level_set/level_set_parameter.h"
#include "level_set/filter_base.h"


MAST::HomogenizedDensityFunctionBase::
HomogenizedDensityFunctionBase(const std::string& nm):
MAST::FieldFunction<Real>  (nm),
_level_set_sys   (nullptr),
_analysis_mesh   (nullptr),
_level_set       (nullptr),
_filter          (nullptr) {
    
}


MAST::HomogenizedDensityFunctionBase::
~HomogenizedDensityFunctionBase() {
    
}


void
MAST::HomogenizedDensityFunctionBase::init(MAST::SystemInitialization& level_set_sys,
                                           libMesh::MeshBase&          analysis_mesh,
                                           MAST::FieldFunction<Real>&  level_set,
                                           MAST::FilterBase&           filter) {
    
    libmesh_assert(!_level_set_sys);
    
    _level_set_sys     = &level_set_sys;
    _analysis_mesh     = &analysis_mesh;
    _level_set         = &level_set;
    _filter            = &filter;
    _sub_point_locator.reset(_analysis_mesh->sub_point_locator().release());
}


void
MAST::HomogenizedDensityFunctionBase::
operator() (const libMesh::Point& p, const Real t, Real& v) const {
    
    libmesh_assert(_analysis_mesh);
    libmesh_assert(!_elem_volume_fraction.empty());
    
    // identify the element in the analysis mesh that this point belongs to
    // and return the value of the homogenized function
    const libMesh::Elem* e = (*_sub_point_locator)(p);

    libmesh_assert(e);
    
    std::map<const libMesh::Elem*, Real>::const_iterator
    it = _elem_volume_fraction.find(e);

    libmesh_assert(it != _elem_volume_fraction.end());

    v = it->second;
}
    


void
MAST::HomogenizedDensityFunctionBase::
derivative(const MAST::FunctionBase& f,
           const libMesh::Point& p, const Real t, Real& v) const {
    
    libmesh_assert(_analysis_mesh);
    libmesh_assert(!_elem_volume_fraction_sensitivity.empty());
    
    // identify the element in the analysis mesh that this point belongs to
    // and return the value of the homogenized function
    const libMesh::Elem* e = (*_sub_point_locator)(p);
    
    libmesh_assert(e);
    
    v = this->get_elem_volume_fraction_sensitivity(f, *e);
}


Real
MAST::HomogenizedDensityFunctionBase::
get_elem_volume_fraction(const libMesh::Elem& e) const {
    
    std::map<const libMesh::Elem*, Real>::const_iterator
    it = _elem_volume_fraction.find(&e);
    
    libmesh_assert(it != _elem_volume_fraction.end());
    
    return it->second;
}



Real
MAST::HomogenizedDensityFunctionBase::
get_elem_volume_fraction_sensitivity(const MAST::FunctionBase& f,
                                     const libMesh::Elem& e) const {
    
    std::map<const MAST::FunctionBase*, std::map<const libMesh::Elem*, Real>>::const_iterator
    it_f   = _elem_volume_fraction_sensitivity.find(&f);

    libmesh_assert (it_f != _elem_volume_fraction_sensitivity.end());

    std::map<const libMesh::Elem*, Real>::const_iterator
    it_e = it_f->second.find(&e);
    
    libmesh_assert (it_e != it_f->second.end());
    
    return it_e->second;
}


const std::map<const libMesh::Elem*, Real>*
MAST::HomogenizedDensityFunctionBase::
get_elem_volume_fraction_sensitivity_map(const MAST::FunctionBase& f) const
{
    std::map<const MAST::FunctionBase*, std::map<const libMesh::Elem*, Real>>::const_iterator
    it   = _elem_volume_fraction_sensitivity.find(&f);
    
    if (it != _elem_volume_fraction_sensitivity.end())
        return &(it->second);
    else
        return nullptr;
}
