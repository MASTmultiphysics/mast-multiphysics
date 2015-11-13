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

#ifndef __mast_beam_optimization_base_h__
#define __mast_beam_optimization_base_h__

// C++ includes
#include <map>

// MAST includes
#include "base/mast_data_types.h"
#include "base/constant_field_function.h"
#include "base/physics_discipline_base.h"


// libMesh includes
#include "libmesh/point.h"

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
        
        MultilinearInterpolation(const MAST::MultilinearInterpolation& o);
        
        
        virtual std::auto_ptr<MAST::FieldFunction<Real> > clone() const;
        
        
        virtual ~MultilinearInterpolation();
        
    protected:
        
        std::map<Real, MAST::FieldFunction<Real>*> _values;
        
    public:
        
        virtual void operator() (const libMesh::Point& p, Real t, Real& v) const;
        
        virtual void derivative(const MAST::DerivativeType d,
                                const MAST::FunctionBase& f,
                                const libMesh::Point& p,
                                Real t,
                                Real& v) const;
    };
    
    
    /*!
     *   Function object evaluates the beam offset for the specified height
     */
    class BeamOffset: public MAST::FieldFunction<Real> {
    public:
        BeamOffset(const std::string& nm,
                   MAST::FieldFunction<Real> *thickness);
        
        BeamOffset(const MAST::BeamOffset& o);
        
        virtual std::auto_ptr<MAST::FieldFunction<Real> > clone() const;
        
        virtual ~BeamOffset();
        
    protected:
        
        MAST::FieldFunction<Real> *_dim;
        
    public:
        
        virtual void operator() (const libMesh::Point& p, Real t, Real& v) const;
        
        virtual void derivative(const MAST::DerivativeType d,
                                const MAST::FunctionBase& f,
                                const libMesh::Point& p,
                                Real t,
                                Real& v) const;
        
    };
    
    
    
    /*!
     *   Function object evaluates the weight and its sensitivity with
     *   respect to the specified variable.
     */
    class Weight: public MAST::FieldFunction<Real> {
    public:
        
        /*!
         *   Constructor requires the mesh and the
         */
        Weight(MAST::PhysicsDisciplineBase& discipline);
        
        /*!
         *  copy constructor
         */
        Weight(const MAST::Weight& w);
        
        /*!
         *  @returns a new object as a clone, encapsulated in a smart-pointer
         */
        virtual std::auto_ptr<MAST::FieldFunction<Real> > clone() const;
        
        /*!
         *  virtual destructor
         */
        virtual ~Weight();
        
    protected:
        
        /*!
         *  discipline object provides the mesh and material properties for
         *  calculation of the mass
         */
        MAST::PhysicsDisciplineBase& _discipline;
        
    public:
        
        /*!
         *   overloaded operator evaluates and returns the mass of the given
         *   structural model.
         */
        virtual void operator() (const libMesh::Point& p,
                                 Real t,
                                 Real& v) const ;
        
        
        
        /*!
         *   evaluates the sensitivity of structural mass with respect
         *   to the design variable.
         */
        virtual void derivative(const MAST::DerivativeType d,
                                const MAST::FunctionBase& f,
                                const libMesh::Point& p,
                                Real t,
                                Real& v) const ;
    };
}


#endif // __mast_beam_optimization_base_h__

