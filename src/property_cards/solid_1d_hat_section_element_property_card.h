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

#ifndef __mast__solid_1d_hat_section_element_property_card__
#define __mast__solid_1d_hat_section_element_property_card__


// MAST includes
#include "property_cards/solid_1d_4parameter_section_element_property_card.h"

// C-Beam C-Channel (CHAN in Nastran)
namespace MAST{
    namespace Solid1DHatSectionProperty{

        void calcA(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& A);

        void calcdA(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dA);

        void calcIz(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& Iz);

        void calcdIz(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dIz);

        void calcIy(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& Iy);

        void calcdIy(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dIy);

        void calcIp(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& Ip);

        void calcdIp(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dIp);

        void calcJ1_w(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& J1_w);

        void calcdJ1_w(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dJ1_w);

        void calcJ2_w(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& J2_w);

        void calcdJ2_w(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dJ2_w);

        void calcJ1_ft(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& J1_ft);

        void calcdJ1_ft(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dJ1_ft);

        void calcJ2_ft(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& J2_ft);

        void calcdJ2_ft(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dJ2_ft);

        void calcJ1_fb(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& J1_fb);

        void calcdJ1_fb(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dJ1_fb);

        void calcJ2_fb(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& J2_fb);

        void calcdJ2_fb(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dJ2_fb);

        void calck1_fb(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& k1_fb);

        void calcdk1_fb(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dk1_fb);

        void calck2_fb(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& k2_fb);

        void calcdk2_fb(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dk2_fb);

        void calck1_ft(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& k1_ft);

        void calcdk1_ft(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dk1_ft);

        void calck2_ft(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& k2_ft);

        void calcdk2_ft(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dk2_ft);

        void calck1_w(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& k1_w);

        void calcdk1_w(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dk1_w);

        void calck2_w(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& k2_w);

        void calcdk2_w(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dk2_w);

        void calcJc(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& Jc);

        void calcdJc(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dJc);
        
        void calcJ(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& J);
        
        void calcdJ(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dJ);

    }
}


namespace MAST {
    
    class Solid1DHatSectionElementPropertyCard : public MAST::Solid1D4ParameterSectionElementPropertyCard
    {
    public:
        Solid1DHatSectionElementPropertyCard():
        MAST::Solid1D4ParameterSectionElementPropertyCard() {}
        
        /*!
         *   virtual destructor
         */
        virtual ~Solid1DHatSectionElementPropertyCard() { }
        
        virtual void init(); 
    };
}


#endif // __mast__solid_1d_hat_section_element_property_card__
