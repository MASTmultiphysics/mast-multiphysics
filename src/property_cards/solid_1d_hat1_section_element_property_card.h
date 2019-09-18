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

#ifndef __mast__solid_1d_hat1_section_element_property_card__
#define __mast__solid_1d_hat1_section_element_property_card__


// MAST includes
#include "property_cards/solid_1d_5parameter_section_element_property_card.h"

// Closed HAT (Hat1 in Nastran and HAT in Astros 21.2)
namespace MAST{
    namespace Solid1DHat1SectionProperty{
        void calcA(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& A);

        void calcdA(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dA);

        void calcIz(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& Iz);

        void calcdIz(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dIz);

        void calcIy(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& Iy);

        void calcdIy(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dIy);

        void calcIp(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& Ip);

        void calcdIp(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dIp);

        void calcJ(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& J);

        void calcdJ(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dJ);

    }
}


namespace MAST {
    
    class Solid1DHat1SectionElementPropertyCard : public MAST::Solid1D5ParameterSectionElementPropertyCard
    {
    public:
        Solid1DHat1SectionElementPropertyCard():
        MAST::Solid1D5ParameterSectionElementPropertyCard() {}
        
        /*!
         *   virtual destructor
         */
        virtual ~Solid1DHat1SectionElementPropertyCard() { }
        
        virtual void init(); 
    };
}


#endif // __mast__solid_1d_hat1_section_element_property_card__
