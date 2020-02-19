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

//TODO: Arbitrary cross sections currently only support sensitivites w.r.t. 
// offsets. Other sensitivies will be returned as zero.

#ifndef __mast__solid_1d_arbitrary_section_element_property_card__
#define __mast__solid_1d_arbitrary_section_element_property_card__


// MAST includes
#include "property_cards/solid_1d_section_element_property_card.h"

namespace MAST {
    
    class Solid1DArbitrarySectionElementPropertyCard : public MAST::Solid1DSectionElementPropertyCard
    {
    public:
        Solid1DArbitrarySectionElementPropertyCard():
        MAST::Solid1DSectionElementPropertyCard() {}
        
        /*!
         *   virtual destructor
         */
        virtual ~Solid1DArbitrarySectionElementPropertyCard() { }
        
        // TODO: Reconsider three below in the source file (.cpp)
        // virtual void init(libMesh::MeshBase& mesh);
        // virtual void init(RealMatrixX& vertices);
        // virtual void init(Real A, Real Izz, Real Iyy, Real J, Real T);
        virtual void init();
        
        void setTorsionalConstant(Real T);
                
        void clear();
        
        virtual const libMesh::Point get_centroid(const libMesh::Point& p, const Real t) const override;
        
        virtual const libMesh::Point get_shear_center(const libMesh::Point& p, const Real t) const override;
        
        virtual const std::vector<libMesh::Point> get_stress_points(const libMesh::Point& p, const Real t, const libMesh::Point ps) const override;
        
        virtual const libMesh::Point get_shear_center_derivative(MAST::FunctionBase& f, const libMesh::Point& p, const Real t) override;
                
        virtual const libMesh::Point get_centroid_derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t) const override;
                
        virtual const std::vector<libMesh::Point> get_stress_points_derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t, const libMesh::Point dps) const override;
        
        virtual void add_stress_point(const libMesh::Point stress_point);
        
    protected:
        
        bool _torsionConstantSet = false;
        
        void calculateGeometricProperties(libMesh::MeshBase& mesh);
        void calculateGeometricProperties(RealMatrixX& vertices);
        
        Real _A_val;
        Real _Izz_val;
        Real _Iyy_val;
        Real _Ip_val;
        Real _J_val;
        Real _W_val;
        Real _Kxx_val;
        Real _Kyy_val;
        std::vector<libMesh::Point> _stress_points;
    };
        
}


#endif // __mast__solid_1d_arbitrary_section_element_property_card__
