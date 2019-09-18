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

#ifndef __mast__solid_1d_pi_section_element_property_card__
#define __mast__solid_1d_pi_section_element_property_card__


// MAST includes
#include "property_cards/solid_1d_5parameter_section_element_property_card.h"

// Pi-Channel (looks like math pi)
namespace MAST{
    namespace Solid1DPiSectionProperty{
        void calcA(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& A);

        void calcdA(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dA);

        void calcIz(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& Iz);

        void calcdIz(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dIz);

        void calcIy(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& Iy);

        void calcdIy(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dIy);

        void calcIp(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& Ip);

        void calcdIp(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dIp);

        void calcJ1_w(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& J1_w);

        void calcdJ1_w(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dJ1_w);

        void calcJ2_w(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& J2_w);

        void calcdJ2_w(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dJ2_w);

        void calcJ1_f(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& J1_f);

        void calcdJ1_f(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dJ1_f);

        void calcJ2_f(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& J2_f);

        void calcdJ2_f(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dJ2_f);

        void calck1_f(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& k1_f);

        void calcdk1_f(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dk1_f);

        void calck2_f(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& k2_f);

        void calcdk2_f(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dk2_f);

        void calck1_w(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& k1_w);

        void calcdk1_w(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dk1_w);

        void calck2_w(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& k2_w);

        void calcdk2_w(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dk2_w);

        void calcJc(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& Jc);

        void calcdJc(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dJc);

        void calcJ(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& J);

        void calcdJ(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dJ);
        
    }
}



namespace MAST {
    
    class Solid1DPiSectionElementPropertyCard : public MAST::Solid1D5ParameterSectionElementPropertyCard
    {
    public:
        Solid1DPiSectionElementPropertyCard():
        MAST::Solid1D5ParameterSectionElementPropertyCard() {}
        
        /*!
         *   virtual destructor
         */
        virtual ~Solid1DPiSectionElementPropertyCard() { }
        
        virtual void init(); 
    };
}


#endif // __mast__solid_1d_pi_section_element_property_card__
