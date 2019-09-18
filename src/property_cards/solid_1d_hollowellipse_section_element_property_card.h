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
//  * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef __mast__solid_1d_hollowellipse_section_element_property_card__
#define __mast__solid_1d_hollowellipse_section_element_property_card__


// MAST includes
#include "property_cards/solid_1d_4parameter_section_element_property_card.h"

// Hollow Rectangle (BOX in Siemens NX Nastran and Astros 21.2)
namespace MAST{
    namespace Solid1DHollowEllipseSectionProperty{
        void calcA(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& A);

        void calcdA(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dA);

        void calcIz(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& Iz);

        void calcdIz(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dIz);

        void calcIy(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& Iy);

        void calcdIy(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dIy);

        void calcIp(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& Ip);

        void calcdIp(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dIp);

        void calcJ(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& J);

        void calcdJ(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dJ);

    }
}


namespace MAST {
    
    class Solid1DHollowEllipseSectionElementPropertyCard : public MAST::Solid1D4ParameterSectionElementPropertyCard
    {
    public:
        Solid1DHollowEllipseSectionElementPropertyCard():
        MAST::Solid1D4ParameterSectionElementPropertyCard() {}
        
        /*!
         *   virtual destructor
         */
        virtual ~Solid1DHollowEllipseSectionElementPropertyCard() { }
        
        virtual void init(); 
    };
}


#endif // __mast__solid_1d_hollowellipse_section_element_property_card__
