///*
// * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
// * Copyright (C) 2013-2020  Manav Bhatia and MAST authors
// *
// * This library is free software; you can redistribute it and/or
// * modify it under the terms of the GNU Lesser General Public
// * License as published by the Free Software Foundation; either
// * version 2.1 of the License, or (at your option) any later version.
// *
// * This library is distributed in the hope that it will be useful,
// * but WITHOUT ANY WARRANTY; without even the implied warranty of
// * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// * Lesser General Public License for more details.
// *
// * You should have received a copy of the GNU Lesser General Public
// * License along with this library; if not, write to the Free Software
// * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
// */
//
//#ifndef __mast_beam_optimization_base_h__
//#define __mast_beam_optimization_base_h__
//
//// C++ includes
//#include <map>
//
//// MAST includes
//#include "base/mast_data_types.h"
//#include "base/constant_field_function.h"
//#include "base/physics_discipline_base.h"
//
//
//// libMesh includes
//#include "libmesh/point.h"
//
//namespace MAST {
//    
//    
//    /*!
//     *   Function object evaluates the BeamWeight and its sensitivity with
//     *   respect to the specified variable.
//     */
//    class BeamWeight: public MAST::FieldFunction<Real> {
//    public:
//        
//        /*!
//         *   Constructor requires the mesh and the
//         */
//        BeamWeight(MAST::PhysicsDisciplineBase& discipline);
//        
//        /*!
//         *  virtual destructor
//         */
//        virtual ~BeamWeight();
//        
//    protected:
//        
//        /*!
//         *  discipline object provides the mesh and material properties for
//         *  calculation of the mass
//         */
//        MAST::PhysicsDisciplineBase& _discipline;
//        
//    public:
//        
//        /*!
//         *   overloaded operator evaluates and returns the mass of the given
//         *   structural model.
//         */
//        virtual void operator() (const libMesh::Point& p,
//                                 Real t,
//                                 Real& v) const ;
//        
//        
//        
//        /*!
//         *   evaluates the sensitivity of structural mass with respect
//         *   to the design variable.
//         */
//        virtual void derivative(   const MAST::FunctionBase& f,
//                                const libMesh::Point& p,
//                                Real t,
//                                Real& v) const ;
//    };
//}
//
//
//#endif // __mast_beam_optimization_base_h__
//
