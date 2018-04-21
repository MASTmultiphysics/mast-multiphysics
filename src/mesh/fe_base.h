/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2018  Manav Bhatia
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


#ifndef __mast_fe_base_h__
#define __mast_fe_base_h__

// MAST includes
#include "base/mast_data_types.h"

// libMesh includes
#include "libmesh/elem.h"
#include "libmesh/fe_base.h"
#include "libmesh/quadrature.h"


namespace MAST {

    // forward declerations
    class SystemInitialization;
    
    class FEBase {
    
    public:
        
        FEBase(const MAST::SystemInitialization& sys);
        
        virtual ~FEBase();
        
        /*!
         *   this is used, in addition to \p  libMesh::System::extra_quadrature_order
         *   to set the quadrature rule.
         */
        void set_extra_quadrature_order(int n);
        
        /*!
         *   Initializes the quadrature and finite element for element volume
         *   integration.
         *   \param elem libMesh::Elem for which the finite element is
         *   initialized. \param pts the points at which the element should
         *   be initialized. If nullptr, the points specified by quadrature
         *   used and the extra quadrature order specified in libMesh::System
         *   is used. An additional change in quadrature order can be specified
         *   using the method set_extra_quadrature_order(). This can be a
         *   positive or negative value.
         */
        virtual void init(const libMesh::Elem& elem,
                          const std::vector<libMesh::Point>* pts = nullptr);

        /*!
         *   Initializes the quadrature and finite element for element side
         *   integration. The last argument
         *   \p if_calculate_dphi tells the function to request the
         *   \p fe object to also initialize the calculation of shape function
         *   derivatives
         */
        virtual void init_for_side(const libMesh::Elem& elem,
                                   unsigned int s,
                                   bool if_calculate_dphi);
        
        libMesh::FEType
        get_fe_type() const;
        
        virtual const std::vector<Real>&
        get_JxW() const;
        
        virtual const std::vector<libMesh::Point>&
        get_xyz() const;

        virtual unsigned int
        n_shape_functions() const;
        
        virtual const std::vector<std::vector<Real> >&
        get_phi() const;

        virtual const std::vector<std::vector<libMesh::RealVectorValue> >&
        get_dphi() const;

        virtual const std::vector<std::vector<libMesh::RealTensorValue>>&
        get_d2phi() const;
        
        virtual const std::vector<Real>&
        get_dxidx() const;

        virtual const std::vector<Real>&
        get_dxidy() const;

        virtual const std::vector<Real>&
        get_dxidz() const;

        virtual const std::vector<Real>&
        get_detadx() const;

        virtual const std::vector<Real>&
        get_detady() const;

        virtual const std::vector<Real>&
        get_detadz() const;

        virtual const std::vector<Real>&
        get_dzetadx() const;

        virtual const std::vector<Real>&
        get_dzetady() const;

        virtual const std::vector<Real>&
        get_dzetadz() const;

        virtual const std::vector<libMesh::RealVectorValue>&
        get_dxyzdxi() const;

        virtual const std::vector<libMesh::RealVectorValue>&
        get_dxyzdeta() const;

        virtual const std::vector<libMesh::RealVectorValue>&
        get_dxyzdzeta() const;

        virtual const std::vector<std::vector<Real> >&
        get_dphidxi() const;

        virtual const std::vector<std::vector<Real> >&
        get_dphideta() const;

        virtual const std::vector<std::vector<Real> >&
        get_dphidzeta() const;

        virtual const std::vector<libMesh::Point>&
        get_normals() const;
        
        virtual const std::vector<libMesh::Point>&
        get_qpoints() const;
        
        virtual const libMesh::QBase&
        get_qrule() const;
        
    protected:
        
        const MAST::SystemInitialization& _sys;
        unsigned int                      _extra_quadrature_order;
        bool                              _init_second_order_derivatives;
        bool                              _initialized;
        const libMesh::Elem*              _elem;
        libMesh::FEBase*                  _fe;
        libMesh::QBase*                   _qrule;
        std::vector<libMesh::Point>       _qpoints;
    };
}


#endif // __mast_fe_base_h__
