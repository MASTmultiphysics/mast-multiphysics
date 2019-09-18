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

#ifndef __mast__solid_1d_chan1_section_element_property_card__
#define __mast__solid_1d_chan1_section_element_property_card__


// MAST includes
#include "property_cards/solid_1d_4parameter_section_element_property_card.h"

// C-Beam C-Channel (CHAN in Nastran)
namespace MAST{
    namespace Solid1DChan1SectionProperty{
        void calcA(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& A);

        void calcdA(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dA);

        void calcIz(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& Iz);

        void calcdIz(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dIz);

        void calcIy(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& Iy);

        void calcdIy(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dIy);

        void calcIp(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& Ip);

        void calcdIp(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dIp);

        void calcJ1_f(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& J1_f);

        void calcdJ1_f(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dJ1_f);

        void calcJ2_f(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& J2_f);

        void calcdJ2_f(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dJ2_f);

        void calcJ1_w(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& J1_w);

        void calcdJ1_w(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dJ1_w);

        void calcJ2_w(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& J2_w);

        void calcdJ2_w(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dJ2_w);

        void calck1_f(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& k1_f);

        void calcdk1_f(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dk1_f);

        void calck2_f(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& k2_f);

        void calcdk2_f(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dk2_f);

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
    
    class Solid1DChan1SectionElementPropertyCard : public MAST::Solid1D4ParameterSectionElementPropertyCard
    {
    public:
        Solid1DChan1SectionElementPropertyCard():
        MAST::Solid1D4ParameterSectionElementPropertyCard() {}
        
        /*!
         *   virtual destructor
         */
        virtual ~Solid1DChan1SectionElementPropertyCard() { }
        
        virtual void init(); 
    };
}


#endif // __mast__solid_1d_chan1_section_element_property_card__
