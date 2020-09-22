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

#ifndef __mast__solid_1d_L_section_element_property_card__
#define __mast__solid_1d_L_section_element_property_card__


// MAST includes
// #include "property_cards/solid_1d_4parameter_section_element_property_card.h"
#include "property_cards/solid_1d_nparameter_section_element_property_card.h"
#include "property_cards/cross_section_property_pilkey.h"


// "L" cross section
namespace MAST{
    namespace Solid1DLSectionProperty{
        void calcA(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& A);

        void calcdA(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dA);

        void calcIz(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& Iz);

        void calcdIz(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dIz);

        void calcIy(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& Iy);

        void calcdIy(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dIy);

        void calcIp(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& Ip);

        void calcdIp(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dIp);
        
        void calcJ1_h(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& J1_h);

        void calcdJ1_h(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dJ1_h);

        void calcJ2_h(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& J2_h);

        void calcdJ2_h(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dJ2_h);

        void calcJ1_v(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& J1_v);

        void calcdJ1_v(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dJ1_v);

        void calcJ2_v(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& J2_v);

        void calcdJ2_v(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dJ2_v);

        void calcJ(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& J);

        void calcdJ(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dJ);

    }
}


namespace MAST {
    
    class Solid1DLSectionElementPropertyCard : public MAST::Solid1DnParameterSectionElementPropertyCard
    {
    public:
        Solid1DLSectionElementPropertyCard():
        MAST::Solid1DnParameterSectionElementPropertyCard() {}
        
        /*!
         *   virtual destructor
         */
        virtual ~Solid1DLSectionElementPropertyCard() { }
        
        virtual void init(const libMesh::LibMeshInit& init,
                          const uint n_target_elems=3500,
                          const libMesh::ElemType element_type=libMesh::TRI6);
        
        virtual void create_cross_section(
            const libMesh::LibMeshInit& init,
            const uint n_target_elems=3500,
            const libMesh::ElemType element_type=libMesh::TRI6);
        
        virtual void calculate_properties_pilkey();
        
        virtual const std::vector<libMesh::Point> get_geom_points(const libMesh::Point& p, const Real t, const uint n=201) const override;
        
        virtual const std::vector<libMesh::Point> get_geom_points_derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t, const uint n=201) const override;
                
        virtual const libMesh::Point get_shear_center(const libMesh::Point& p, const Real t) const override;
        
        virtual const libMesh::Point get_shear_center_derivative(MAST::FunctionBase& f, const libMesh::Point& p, const Real t) override;
        
        virtual const libMesh::Point get_centroid(const libMesh::Point& p, const Real t) const override;
        
        virtual const libMesh::Point get_centroid_derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t) const override;
        
        virtual const std::vector<libMesh::Point> get_stress_points(const libMesh::Point& p, const Real t, const libMesh::Point ps) const override;
        
        virtual const std::vector<libMesh::Point> get_stress_points_derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t, const libMesh::Point dps) const override;
        
        std::unique_ptr<MAST::CrossSection> cross_section;
        
        const bool bysymmetric = false;
        const bool symmetric_z = false;
        const bool symmetric_y = false;
    };
}


#endif // __mast__solid_1d_L_section_element_property_card__
