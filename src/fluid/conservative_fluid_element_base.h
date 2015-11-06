/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2015  Manav Bhatia
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

#ifndef __mast__conservative_fluid_element_base__
#define __mast__conservative_fluid_element_base__


// MAST includes
#include "base/elem_base.h"
#include "fluid/fluid_elem_base.h"


namespace MAST {
    
    // Forward declerations
    class ElementPropertyCardBase;
    class LocalElemBase;
    class BoundaryConditionBase;
    class FEMOperatorMatrix;

    
    /*!
     *  This class provides the necessary functionality for spatial
     *  discretization of the conservative fluid equations.
     */
    class ConservativeFluidElementBase:
    public MAST::FluidElemBase,
    public MAST::ElementBase {
        
    public:
        
        ConservativeFluidElementBase(MAST::SystemInitialization& sys,
                                     const libMesh::Elem& elem,
                                     const MAST::ElementPropertyCardBase& p);
        
        virtual ~ConservativeFluidElementBase();
        
        /*!
         *   returns a constant reference to the finite element object
         */
        const MAST::ElementPropertyCardBase&
        elem_property()  {
            return _property;
        }
        
        
        /*!
         *   internal force contribution to system residual
         */
        virtual bool
        internal_residual (bool request_jacobian,
                           RealVectorX& f,
                           RealMatrixX& jac);
        
        
        /*!
         *   inertial force contribution to system residual
         */
        virtual bool
        velocity_residual (bool request_jacobian,
                           RealVectorX& f,
                           RealMatrixX& jac_xdot,
                           RealMatrixX& jac);
        
        /*!
         *   side external force contribution to system residual
         */
        bool
        side_external_residual (bool request_jacobian,
                                RealVectorX& f,
                                RealMatrixX& jac,
                                std::multimap<libMesh::boundary_id_type, MAST::BoundaryConditionBase*>& bc);
        
        /*!
         *   sensitivity of the internal force contribution to system residual
         */
        virtual bool
        internal_residual_sensitivity (bool request_jacobian,
                                       RealVectorX& f,
                                       RealMatrixX& jac);
        /*!
         *   sensitivity of the damping force contribution to system residual
         */
        virtual bool
        velocity_residual_sensitivity (bool request_jacobian,
                                       RealVectorX& f,
                                       RealMatrixX& jac);
        
        /*!
         *   sensitivity of the side external force contribution to system residual
         */
        bool
        side_external_residual_sensitivity (bool request_jacobian,
                                            RealVectorX& f,
                                            RealMatrixX& jac,
                                            std::multimap<libMesh::boundary_id_type, MAST::BoundaryConditionBase*>& bc);
        
        /*!
         *   evaluates an output quantity requested in the map over the
         *   boundary of the element that may coincide with the boundary
         *   identified in the map. The derivative with respect to the
         *   state variables is provided if \p request_derivative is true.
         */
        virtual bool
        volume_output_quantity (bool request_derivative,
                                bool request_sensitivity,
                                std::multimap<libMesh::subdomain_id_type, MAST::OutputFunctionBase*>& output) {
            
            libmesh_error(); // to be implemented
        }
        
        
        /*!
         *   evaluates an output quantity requested in the map over the
         *   boundary of the element that may coincide with the boundary
         *   identified in the map. The derivative with respect to the
         *   state variables is provided if \p request_derivative is true.
         */
        virtual bool
        side_output_quantity (bool request_derivative,
                              std::multimap<libMesh::boundary_id_type, MAST::OutputFunctionBase*>& output) {
            
            libmesh_error(); // to be implemented
        }

        
    protected:
        
        
        
        /*!
         *
         */
        virtual bool symmetry_surface_residual(bool request_jacobian,
                                               RealVectorX& f,
                                               RealMatrixX& jac,
                                               const unsigned int s,
                                               MAST::BoundaryConditionBase& p);
        
        /*!
         *
         */
        virtual bool symmetry_surface_residual_sensitivity(bool request_jacobian,
                                                           RealVectorX& f,
                                                           RealMatrixX& jac,
                                                           const unsigned int s,
                                                           MAST::BoundaryConditionBase& p);
        
        /*!
         *
         */
        virtual bool slip_wall_surface_residual(bool request_jacobian,
                                                RealVectorX& f,
                                                RealMatrixX& jac,
                                                const unsigned int s,
                                                MAST::BoundaryConditionBase& p);
        
        /*!
         *
         */
        virtual bool slip_wall_surface_residual_sensitivity(bool request_jacobian,
                                                            RealVectorX& f,
                                                            RealMatrixX& jac,
                                                            const unsigned int s,
                                                            MAST::BoundaryConditionBase& p);
        
        /*!
         *
         */
        virtual bool far_field_surface_residual(bool request_jacobian,
                                                RealVectorX& f,
                                                RealMatrixX& jac,
                                                const unsigned int s,
                                                MAST::BoundaryConditionBase& p);
        
        /*!
         *
         */
        virtual bool far_field_surface_residual_sensitivity(bool request_jacobian,
                                                            RealVectorX& f,
                                                            RealMatrixX& jac,
                                                            const unsigned int s,
                                                            MAST::BoundaryConditionBase& p);
        

        /*!
         */
        void _initialize_fem_interpolation_operator(const unsigned int qp,
                                                    const unsigned int dim,
                                                    const libMesh::FEBase& fe,
                                                    MAST::FEMOperatorMatrix& Bmat);
        
        
        /*!
         *    For \p mass = true, the FEM operator matrix is initilized to
         *    the weak form of the Laplacian
         *    \f[  dB[0] = \frac{\partial {\bf N}}{\partial x} \f]
         *
         *    \f[  dB[1] = \frac{\partial {\bf N}}{\partial y} \f]
         *
         *    \f[  dB[2] = \frac{\partial {\bf N}}{\partial z} \f]
         */
        void _initialize_fem_gradient_operator(const unsigned int qp,
                                               const unsigned int dim,
                                               const libMesh::FEBase& fe,
                                               std::vector<MAST::FEMOperatorMatrix>& dBmat);
        
        /*!
         *   element property
         */
        const MAST::ElementPropertyCardBase& _property;
        
    };
}


#endif // __mast__conservative_fluid_element_base__
