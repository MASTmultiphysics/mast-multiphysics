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

#ifndef __mast__solid_1d_I_section_element_property_card__
#define __mast__solid_1d_I_section_element_property_card__


// MAST includes
#include "property_cards/solid_1d_6parameter_section_element_property_card.h"

// I-Beam Asymmetric about z-axis (I in Nastran and HAT in Astros 21.2)
namespace MAST{
    namespace Solid1DISectionProperty{
        void calcA(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& A);

        void calcdA(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dDIM6, Real& dA);

        void calcIz(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& Iz);

        void calcdIz(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dDIM6, Real& dIz);

        void calcIy(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& Iy);

        void calcdIy(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dDIM6, Real& dIy);

        void calcIp(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& Ip);

        void calcdIp(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dDIM6, Real& dIp);

        void calcJ1_w(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& J1_w);

        void calcdJ1_w(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dDIM6, Real& dJ1_w);

        void calcJ2_w(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& J2_w);

        void calcdJ2_w(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dDIM6, Real& dJ2_w);

        void calcJ1_ft(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& J1_ft);

        void calcdJ1_ft(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dDIM6, Real& dJ1_ft);

        void calcJ2_ft(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& J2_ft);

        void calcdJ2_ft(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dDIM6, Real& dJ2_ft);

        void calcJ1_fb(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& J1_fb);

        void calcdJ1_fb(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dDIM6, Real& dJ1_fb);

        void calcJ2_fb(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& J2_fb);

        void calcdJ2_fb(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dDIM6, Real& dJ2_fb);

        void calck1_ft(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& k1_ft);

        void calcdk1_ft(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dDIM6, Real& dk1_ft);

        void calck2_ft(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& k2_ft);

        void calcdk2_ft(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dDIM6, Real& dk2_ft);

        void calck1_fb(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& k1_fb);

        void calcdk1_fb(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dDIM6, Real& dk1_fb);

        void calck2_fb(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& k2_fb);

        void calcdk2_fb(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dDIM6, Real& dk2_fb);

        void calck1_w(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& k1_w);

        void calcdk1_w(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dDIM6, Real& dk1_w);

        void calck2_w(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& k2_w);

        void calcdk2_w(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dDIM6, Real& dk2_w);

        void calcJct(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& Jct);

        void calcdJct(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dDIM6, Real& dJct);

        void calcJcb(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& Jcb);

        void calcdJcb(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dDIM6, Real& dJcb);
        
        void calcJ(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& Jct);

        void calcdJ(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dDIM6, Real& dJct);
    }
}
        




namespace MAST {
    
    class Solid1DISectionElementPropertyCard : public MAST::Solid1D6ParameterSectionElementPropertyCard
    {
    public:
        Solid1DISectionElementPropertyCard():
        MAST::Solid1D6ParameterSectionElementPropertyCard() {}
        
        /*!
         *   virtual destructor
         */
        virtual ~Solid1DISectionElementPropertyCard() { }
        
        virtual void init(); 
    };
}


#endif // __mast__solid_1d_I_section_element_property_card__
