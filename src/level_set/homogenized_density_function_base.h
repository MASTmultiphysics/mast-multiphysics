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

#ifndef __mast__homogenized_density_function_base_h__
#define __mast__homogenized_density_function_base_h__

// MAST includes
#include "base/field_function_base.h"

// libMesh includes
#include "libmesh/mesh_base.h"


namespace MAST {
 
    // Forward declerations
    class FilterBase;
    class SystemInitialization;
    
    
    class HomogenizedDensityFunctionBase:
    public MAST::FieldFunction<Real> {
      
    public:
        
        HomogenizedDensityFunctionBase(const std::string& nm);
        
        virtual ~HomogenizedDensityFunctionBase();

        virtual void init(MAST::SystemInitialization& level_set_sys,
                          libMesh::MeshBase&          analysis_mesh,
                          MAST::FieldFunction<Real>&  level_set,
                          MAST::FilterBase&           filter);
        
        virtual bool depends_on(const MAST::FunctionBase& f) const { return true;}
        
        virtual void operator() (const libMesh::Point& p, const Real t, Real& v) const;
        
        virtual void derivative(const MAST::FunctionBase& f,
                                const libMesh::Point& p, const Real t, Real& v) const;

        const std::map<const libMesh::Elem*, Real>&
        get_elem_volume_fraction_map() const { return _elem_volume_fraction;}
        
        const std::map<const libMesh::Elem*, Real>*
        get_elem_volume_fraction_sensitivity_map(const MAST::FunctionBase& f) const;

        Real
        get_elem_volume_fraction(const libMesh::Elem& e) const;

        Real
        get_elem_volume_fraction_sensitivity(const MAST::FunctionBase& f,
                                             const libMesh::Elem& e) const;

        virtual void initialize_element_volume_fractions() = 0;
        
        virtual void
        initialize_element_volume_fraction_sensitivity(const MAST::FunctionBase& f) = 0;
        
        virtual void
        clear_element_volume_fractions() { _elem_volume_fraction.clear(); }

        virtual void
        clear_element_volume_fraction_sensitivity() { _elem_volume_fraction_sensitivity.clear();}
        
    protected:

        MAST::SystemInitialization           *_level_set_sys;
        
        libMesh::MeshBase                    *_analysis_mesh;
        
        MAST::FieldFunction<Real>            *_level_set;
        
        MAST::FilterBase                     *_filter;
        
        std::map<const libMesh::Elem*, Real> _elem_volume_fraction;
        
        std::map<const MAST::FunctionBase*, std::map<const libMesh::Elem*, Real>>
        _elem_volume_fraction_sensitivity;

        std::unique_ptr<libMesh::PointLocatorBase> _sub_point_locator;
    };
}

#endif // __mast__homogenized_density_function_base_h__

