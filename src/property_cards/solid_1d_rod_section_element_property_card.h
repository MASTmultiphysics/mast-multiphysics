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

#ifndef __mast__solid_1d_rod_section_element_property_card__
#define __mast__solid_1d_rod_section_element_property_card__


// MAST includes
#include "property_cards/solid_1d_1parameter_section_element_property_card.h"


namespace MAST 
{
    namespace Solid1DRodSectionProperty 
    {
        void calcA(Real& DIM1, Real& A);
        
        void calcdA(Real& DIM1, Real& dDIM1, Real& dA);
        
        void calcIz(Real& DIM1, Real& Iz);
        
        void calcdIz(Real& DIM1, Real& dDIM1, Real& dIz);
        
        void calcIy(Real& DIM1, Real& Iy);
        
        void calcdIy(Real& DIM1, Real& dDIM1, Real& dIy);
        
        void calcIp(Real& DIM1, Real& Ip);
        
        void calcdIp(Real& DIM1, Real& dDIM1, Real& dIp);
        
        void calcJ(Real& DIM1, Real& J);
        
        void calcdJ(Real& DIM1, Real& dDIM1, Real& dJ);
        
        class WarpingConstant: public MAST::FieldFunction<Real> 
        {
        public:
            WarpingConstant();
            
            virtual ~WarpingConstant() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const;
            
            virtual void derivative (const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const;
                                     
            /*!
             * Added to be compatible with other sections where the derivative
             * function is not constant due to using finite differences to 
             * approximate the sensitivity.
             */
            virtual void derivative (MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     Real& m) ;
        };
    }
}

namespace MAST {
    
    class Solid1DRodSectionElementPropertyCard : public MAST::Solid1DnParameterSectionElementPropertyCard
    {
    public:
        Solid1DRodSectionElementPropertyCard():
        MAST::Solid1DnParameterSectionElementPropertyCard() {}
        
        /*!
         *   virtual destructor
         */
        virtual ~Solid1DRodSectionElementPropertyCard() { }
        
        virtual void init(const libMesh::LibMeshInit& init); 
        
        virtual void create_cross_section(
            const libMesh::LibMeshInit& init,
            const uint n_target_elems=3500,
            const libMesh::ElemType element_type=libMesh::TRI6);
        
        virtual void calculate_properties_pilkey();
        
        virtual const std::vector<libMesh::Point> get_geom_points(const libMesh::Point& p, const Real t, const uint n=201) const override;
        
        virtual const std::vector<libMesh::Point> get_geom_points_derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t, const uint n=201) const override;
                
        virtual const libMesh::Point get_shear_center(const libMesh::Point& p, const Real t) const override;
        
        virtual const libMesh::Point get_shear_center_derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t) const;
        
        virtual const libMesh::Point get_centroid(const libMesh::Point& p, const Real t) const override;
        
        virtual const libMesh::Point get_centroid_derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t) const override;
        
        virtual const std::vector<libMesh::Point> get_stress_points(const libMesh::Point& p, const Real t, const libMesh::Point ps) const override;
        
        virtual const std::vector<libMesh::Point> get_stress_points_derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t, const libMesh::Point dps) const override;
        
        std::unique_ptr<MAST::CrossSection> cross_section;
        
        const bool is_bisymmetric() const
        {
            return true;
        }
        
        const bool is_symmetric_z() const
        {
            return true;
        }
        
        const bool is_symmetric_y() const
        {
            return true;
        }
    };
}


#endif // __mast__solid_1d_rod_section_element_property_card__
