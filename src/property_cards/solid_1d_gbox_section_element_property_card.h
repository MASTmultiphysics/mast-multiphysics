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

#ifndef __mast__solid_1d_gbox_section_element_property_card__
#define __mast__solid_1d_gbox_section_element_property_card__


// MAST includes
#include "property_cards/solid_1d_6parameter_section_element_property_card.h"


// GBOX Section (GBOX in Astros 21.2, looks like Roman Numeral II)
namespace MAST{
    namespace Solid1DGBOXSectionProperty{
        void calcA(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& A);

        void calcdA(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dDIM6, Real& dA);

        void calcIz(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& Iz);

        void calcdIz(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dDIM6, Real& dIz);

        void calcIy(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& Iy);

        void calcdIy(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dDIM6, Real& dIy);

        void calcIp(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& Ip);

        void calcdIp(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dDIM6, Real& dIp);

        void calcJ_box(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& J_box);

        void calcdJ_box(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dDIM6, Real& dJ_box);

        void calcJ1_fendst(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& J1_fendst);

        void calcdJ1_fendst(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dDIM6, Real& dJ1_fendst);

        void calcJ2_fendst(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& J2_fendst);

        void calcdJ2_fendst(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dDIM6, Real& dJ2_fendst);

        void calcJ1_fendsb(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& J1_fendsb);

        void calcdJ1_fendsb(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dDIM6, Real& dJ1_fendsb);

        void calcJ2_fendsb(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& J2_fendsb);

        void calcdJ2_fendsb(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dDIM6, Real& dJ2_fendsb);

        void calcJc(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& Jc);

        void calcdJc(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dDIM6, Real& dJc);
        
        void calcJ(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& Jc);

        void calcdJ(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dDIM6, Real& dJc);

    }
}

     
namespace MAST {
    
    class Solid1DGboxSectionElementPropertyCard : public MAST::Solid1D6ParameterSectionElementPropertyCard
    {
    public:
        Solid1DGboxSectionElementPropertyCard():
        MAST::Solid1D6ParameterSectionElementPropertyCard() {}
        
        /*!
         *   virtual destructor
         */
        virtual ~Solid1DGboxSectionElementPropertyCard() { }
        
        virtual void init(); 
    };
}


#endif // __mast__solid_1d_gbox_section_element_property_card__
