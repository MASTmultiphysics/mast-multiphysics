/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2019  Manav Bhatia
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

#ifndef __mast_multilinear_interpolation_h__
#define __mast_multilinear_interpolation_h__


// C++ includes
#include <map>

// MAST includes
#include "base/field_function_base.h"



namespace MAST {
    
    
    
    /*!
     *   This class provides the ability to interpolate a function in between
     *   a set of tabulated points.
     */
    class MultilinearInterpolation:
    public MAST::FieldFunction<Real> {
    public:
        MultilinearInterpolation(const std::string& nm,
                                 std::map<Real, MAST::FieldFunction<Real>*>& values);
        
        
        virtual ~MultilinearInterpolation();
        
    protected:
        
        std::map<Real, MAST::FieldFunction<Real>*> _values;
        
    public:
        
        virtual void operator() (const libMesh::Point& p, Real t, Real& v) const;
        
        virtual void derivative(   const MAST::FunctionBase& f,
                                const libMesh::Point& p,
                                Real t,
                                Real& v) const;
    };
    
    
    
    
    /*!
     *   Function object evaluates the beam offset for the specified height
     */
    class SectionOffset: public MAST::FieldFunction<Real> {
    public:
        SectionOffset(const std::string& nm,
                      const MAST::FieldFunction<Real> &thickness,
                      const Real scale);
        
        virtual ~SectionOffset();
        
    protected:
        
        const MAST::FieldFunction<Real> &_dim;
        
        Real _scale;
        
    public:
        
        virtual void operator() (const libMesh::Point& p, Real t, Real& v) const;
        
        virtual void derivative(   const MAST::FunctionBase& f,
                                const libMesh::Point& p,
                                Real t,
                                Real& v) const;
        
    };
    
}

#endif // __mast_multilinear_interpolation_h__

