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

//TODO: Arbitrary cross sections currently only support sensitivites w.r.t. 
// offsets. Other sensitivies will be returned as zero.

#ifndef __mast__solid_1d__arbitrary_section_element_property_card__
#define __mast__solid_1d_arbitrary_section_element_property_card__


// MAST includes
#include "property_cards/solid_1d_section_element_property_card.h"

namespace MAST {
    
    class Solid1DArbitrarySectionElementPropertyCard : public MAST::Solid1DSectionElementPropertyCard
    {
    public:
        Solid1DArbitrarySectionElementPropertyCard():
        MAST::Solid1DSectionElementPropertyCard() {}
        
        /*!
         *   virtual destructor
         */
        virtual ~Solid1DArbitrarySectionElementPropertyCard() { }
        
        virtual void init(libMesh::MeshBase& mesh);
        virtual void init(RealMatrixX& vertices);
        virtual void init(Real A, Real Izz, Real Iyy, Real J, Real T);
        
        void setTorsionalConstant(Real T);
                
        void clear();
        
    protected:
        
        bool _torsionConstantSet = false;
        
        void calculateGeometricProperties(libMesh::MeshBase& mesh);
        void calculateGeometricProperties(RealMatrixX& vertices);
        
        Real _A_val;
        Real _Izz_val;
        Real _Iyy_val;
        Real _Ip_val;
        Real _J_val;
    };
        
}


#endif // __mast__solid_1d_arbitrary_section_element_property_card__
