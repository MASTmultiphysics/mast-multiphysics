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

#ifndef __mast_sub_cell_fe_h__
#define __mast_sub_cell_fe_h__

// MAST includes
#include "mesh/fe_base.h"


namespace MAST {
    
    // Forward declerations
    class LevelSetIntersection;
    
    /*!
     *   This class specializes the MAST::FEBase class for level-set applications
     *   where integration is to be performed on a sub-cell inside an element.
     *   This requires that the quadrature be initialized using the points
     *   in the sub-cell, but evaluated on the parent element. This way, the
     *   shape functions of the original element are used, and the quadrature
     *   weights (JxW) are provided by the sub-element.
     */
    class SubCellFE:
    public MAST::FEBase {
      
    public:
        
        SubCellFE(const MAST::SystemInitialization& sys,
                  const MAST::LevelSetIntersection& intersection);
        
        virtual ~SubCellFE();

        /*!
         *   This assumes that \p elem is the sub-cell and that the original
         *   element will be available as the parent of this element.
         */
        virtual void init(const MAST::GeomElem& elem,
                          bool init_grads,
                          const std::vector<libMesh::Point>* pts = nullptr);
        
        /*!
         *   This assumes that \p elem is the sub-cell and that the original
         *   element will be available as the parent of this element.
         */
        virtual void init_for_side(const MAST::GeomElem& elem,
                                   unsigned int s,
                                   bool if_calculate_dphi);

        virtual const std::vector<Real>&
        get_JxW() const;

        virtual const std::vector<libMesh::Point>&
        get_qpoints() const;

    protected:

      
        const MAST::LevelSetIntersection& _intersection;

        std::vector<Real>   _subcell_JxW;
    };
}

#endif // __mast_sub_cell_fe_h__
