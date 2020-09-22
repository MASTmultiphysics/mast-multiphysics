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


#ifndef __mast__solid_1d_nparameter_section_element_property_card__
#define __mast__solid_1d_nparameter_section_element_property_card__


// MAST includes
#include "property_cards/solid_1d_section_element_property_card.h"
#include "base/field_function_base.h"
#include "property_cards/cross_section_property_pilkey.h"


namespace MAST {
    namespace Solid1DnParameterSectionProperty {
        
        class Area: public MAST::FieldFunction<Real> 
        {
        public:
            Area(MAST::CrossSection& cross_section);
            
            virtual ~Area() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const;
            
            virtual void derivative (const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const;
            
        protected:
            MAST::CrossSection&  _cross_section;
        };
        
        
        /*!
         *   calculates the area moment about the Y-axis due to an offset 
         *   along the Z-axis
         */
        class AreaYMoment: public MAST::FieldFunction<Real> 
        {
        public:
            AreaYMoment(MAST::CrossSection& cross_section);
            
            
            virtual ~AreaYMoment() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const;
            
            
            virtual void derivative (const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const;
            
        protected:
            MAST::CrossSection&  _cross_section;
        };

        
        /*!
         *   calculates the area moment about the Z-axis due to an offset
         *   along the Y-axis
         */
        class AreaZMoment: public MAST::FieldFunction<Real> 
        {
        public:
            AreaZMoment(MAST::CrossSection& cross_section);
            
            virtual ~AreaZMoment() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const;
            
            
            virtual void derivative (const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const;
            
        protected:
            MAST::CrossSection&  _cross_section;
        };
        
        
        /*!
         *   calculates the 2x2 matrix of area inertia for the section
         */
        class AreaInertiaMatrix: public MAST::FieldFunction<RealMatrixX> 
        {
        public:
            AreaInertiaMatrix(MAST::CrossSection& cross_section);
            
            virtual ~AreaInertiaMatrix() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
            
            virtual void derivative (const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
        protected:
            MAST::CrossSection&  _cross_section;
        };
        
        
        class PolarInertia: public MAST::FieldFunction<Real> 
        {
        public:
            PolarInertia(MAST::CrossSection& cross_section);
            
            virtual ~PolarInertia() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const;
            
            virtual void derivative (const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const;
            
        protected:
            MAST::CrossSection&  _cross_section;
        };
        
    
        class TorsionalConstant: public MAST::FieldFunction<Real> 
        {
        public:
            TorsionalConstant(MAST::CrossSection& cross_section);
            
            virtual ~TorsionalConstant() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const;
            
            virtual void derivative (MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     Real& m);
            
        protected:
            MAST::CrossSection&  _cross_section;
        };
        
        
        class WarpingConstant: public MAST::FieldFunction<Real> 
        {
        public:
            WarpingConstant(MAST::CrossSection& cross_section);
            
            virtual ~WarpingConstant() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const;
            
            virtual void derivative (MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     Real& m);
            
        protected:
            MAST::CrossSection&  _cross_section;
        };
        
        
        class ShearCoefficientMatrix: public MAST::FieldFunction<RealMatrixX> 
        {
        public:
            ShearCoefficientMatrix(MAST::CrossSection& cross_section);
            
            virtual ~ShearCoefficientMatrix() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
            
            virtual void derivative (MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m);
            
        protected:
            MAST::CrossSection&  _cross_section;
        };
        
    }
}


namespace MAST {
    
    class Solid1DnParameterSectionElementPropertyCard : public MAST::Solid1DSectionElementPropertyCard
    {
    public:
        Solid1DnParameterSectionElementPropertyCard():
        MAST::Solid1DSectionElementPropertyCard() {}
        
        /*!
         *   virtual destructor
         */
        virtual ~Solid1DnParameterSectionElementPropertyCard() { }
    };
}


#endif // __mast__solid_1d_nparameter_section_element_property_card__
