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

#ifndef __mast_stiffened_plate_optimization_base_h__
#define __mast_stiffened_plate_optimization_base_h__

// C++ includes
#include <map>

// MAST includes
#include "base/mast_data_types.h"
#include "base/constant_field_function.h"
#include "base/physics_discipline_base.h"


// libMesh includes
#include "libmesh/point.h"
#include "libmesh/mesh_base.h"

namespace MAST {
    
    
    /*!
     *   Function object evaluates the PlateWeight and its sensitivity with
     *   respect to the specified variable.
     */
    class StiffenedPlateWeight: public MAST::FieldFunction<Real> {
    public:
        
        /*!
         *   Constructor requires the mesh and the
         */
        StiffenedPlateWeight(MAST::PhysicsDisciplineBase& discipline);
        
        /*!
         *  copy constructor
         */
        StiffenedPlateWeight(const MAST::StiffenedPlateWeight& w);
        
        /*!
         *  @returns a new object as a clone, encapsulated in a smart-pointer
         */
        virtual std::auto_ptr<MAST::FieldFunction<Real> > clone() const;
        
        /*!
         *  virtual destructor
         */
        virtual ~StiffenedPlateWeight();
        
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
    

    /*!
     *    builds the mesh for a stiffened panel
     */
    class StiffenedPanelMesh {
    public:
        StiffenedPanelMesh() { }
        
        ~StiffenedPanelMesh() { }
        
        
        void init (const unsigned int n_stiff,
                   const unsigned int n_x_divs,
                   const unsigned int n_y_divs_between_stiffeners,
                   const Real length,
                   const Real width,
                   libMesh::MeshBase& mesh,
                   libMesh::ElemType t,
                   bool beam_stiffeners);
        
    protected:
        
        enum Component {
            PANEL,
            STIFFENER_X,
            STIFFENER_Y
        };
        
        void _combine_mesh(libMesh::MeshBase& panel,
                           libMesh::MeshBase& stiffener,
                           MAST::StiffenedPanelMesh::Component c,
                           Real stiff_offset,
                           libMesh::subdomain_id_type sid);
        
    };
    
}


#endif // __mast_stiffened__plate_optimization_base_h__

