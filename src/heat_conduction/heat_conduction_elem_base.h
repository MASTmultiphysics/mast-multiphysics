/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2016  Manav Bhatia
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


#ifndef __mast__heat_conduction_elem_base__
#define __mast__heat_conduction_elem_base__

// MAST includes
#include "base/elem_base.h"

namespace MAST {
    
    // Forward declerations
    class ElementPropertyCardBase;
    class LocalElemBase;
    class BoundaryConditionBase;
    class FEMOperatorMatrix;
    
    
    /*!
     *    This element implements the Galerkin discretization of the 
     *    heat conduction problem
     *    \f[ \rho c_p \frac{\partial T}{\partial t} - \frac{\partial }{\partial x_i}\left( -k_{ij} \frac{\partial T}{\partial x_j} \right) = q_v \mbox{ in } \Omega \f]
     *    with the flux provided on the boundary with Neumann boundary conditions. 
     *    The discrete form is represented as
     *    \f[ M^k \dot{x}^k_t + f(x^k_t) = 0 \f],
     *    where 
     *    \f[ M^k \dot{x}^k_t  = \int_{\Omega_e} \phi \rho c_p \dot{x}^k_t \f],
     *    and 
     *    \f[ f(x^k_t) = 
     *    \int_{\Omega_e} \frac{\partial \phi}{\partial x_i} k_ij \frac{\partial T}{\partial x_j}
     *    - \int_{\partial\Omega_e \cap \partial_N\Omega} \phi q_n 
     *    - \int_{\Omega_e} \phi q_v \f],
     *    where \f$q_n\f$ is the boundary normal surface flux on
     *    \f$\partial_N\Omega\f$ and \f$q_v\f$ is the volumetric heat generation.
     */
    class HeatConductionElementBase:
    public MAST::ElementBase
    {
    public:
        /*!
         *   Constructor
         */
        HeatConductionElementBase(MAST::SystemInitialization& sys,
                                  const libMesh::Elem& elem,
                                  const MAST::ElementPropertyCardBase& p);
        
        
        virtual ~HeatConductionElementBase();
        
        
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
         *   volume external force contribution to system residual
         */
        bool
        volume_external_residual (bool request_jacobian,
                                  RealVectorX& f,
                                  RealMatrixX& jac,
                                  std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*>& bc);
        
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
         *   sensitivity of the volume external force contribution to system residual
         */
        bool
        volume_external_residual_sensitivity (bool request_jacobian,
                                              RealVectorX& f,
                                              RealMatrixX& jac,
                                              std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*>& bc);
        
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
         *    Calculates the residual vector and Jacobian due to surface flux
         *    on element side \par s.
         */
        virtual bool surface_flux_residual(bool request_jacobian,
                                           RealVectorX& f,
                                           RealMatrixX& jac,
                                           const unsigned int s,
                                           MAST::BoundaryConditionBase& p);
        
        /*!
         *    Calculates the residual vector sensitivity due to
         *    surface flux on element volumetric domain. This is used only for
         *    2D or 3D elements.
         */
        virtual bool surface_flux_residual(bool request_jacobian,
                                           RealVectorX& f,
                                           RealMatrixX& jac,
                                           MAST::BoundaryConditionBase& p);
        
        /*!
         *    Calculates the residual vector sensitivity and Jacobian due to
         *    surface flux on element side \par s.
         */
        virtual bool surface_flux_residual_sensitivity(bool request_jacobian,
                                                       RealVectorX& f,
                                                       RealMatrixX& jac,
                                                       const unsigned int s,
                                                       MAST::BoundaryConditionBase& p);
        
        /*!
         *    Calculates the residual vector and Jacobian due to surface flux
         *    on element volumetric domain. This is used only for 1D or 2D
         *    elements.
         */
        virtual bool surface_flux_residual_sensitivity(bool request_jacobian,
                                                       RealVectorX& f,
                                                       RealMatrixX& jac,
                                                       MAST::BoundaryConditionBase& p);
        
        
        /*!
         *    Calculates the residual vector and Jacobian due to surface
         *    convection.
         */
        virtual bool surface_convection_residual(bool request_jacobian,
                                                 RealVectorX& f,
                                                 RealMatrixX& jac,
                                                 const unsigned int s,
                                                 MAST::BoundaryConditionBase& p);

        
        /*!
         *    Calculates the residual vector and Jacobian due to surface
         *    convection on the element domain. This is relevant only for
         *    1D and 2D elements.
         */
        virtual bool surface_convection_residual(bool request_jacobian,
                                                 RealVectorX& f,
                                                 RealMatrixX& jac,
                                                 MAST::BoundaryConditionBase& p);

        
        /*!
         *    Calculates the residual vector sensitivity and Jacobian due
         *    to surface convection.
         */
        virtual bool
        surface_convection_residual_sensitivity(bool request_jacobian,
                                                RealVectorX& f,
                                                RealMatrixX& jac,
                                                const unsigned int s,
                                                MAST::BoundaryConditionBase& p);

        /*!
         *    Calculates the residual vector sensitivity and Jacobian due
         *    to surface convection on element domain. This is relevant only 
         *    for 1D and 2D elements.
         */
        virtual bool
        surface_convection_residual_sensitivity(bool request_jacobian,
                                                RealVectorX& f,
                                                RealMatrixX& jac,
                                                MAST::BoundaryConditionBase& p);

        
        /*!
         *    Calculates the residual vector and Jacobian due to surface
         *    radiation flux on side s.
         */
        virtual bool surface_radiation_residual(bool request_jacobian,
                                                RealVectorX& f,
                                                RealMatrixX& jac,
                                                const unsigned int s,
                                                MAST::BoundaryConditionBase& p);

        
        /*!
         *    Calculates the residual vector and Jacobian due to surface 
         *    radiation flux on element domain.
         */
        virtual bool surface_radiation_residual(bool request_jacobian,
                                                RealVectorX& f,
                                                RealMatrixX& jac,
                                                MAST::BoundaryConditionBase& p);

        
        /*!
         *    Calculates the residual vector sensitivity and Jacobian due
         *    to surface radiation flux on element side.
         */
        virtual bool
        surface_radiation_residual_sensitivity(bool request_jacobian,
                                               RealVectorX& f,
                                               RealMatrixX& jac,
                                               const unsigned int s,
                                               MAST::BoundaryConditionBase& p);

        
        /*!
         *    Calculates the residual vector sensitivity and Jacobian due
         *    to surface radiation flux on element domain.
         */
        virtual bool
        surface_radiation_residual_sensitivity(bool request_jacobian,
                                               RealVectorX& f,
                                               RealMatrixX& jac,
                                               MAST::BoundaryConditionBase& p);

        
        /*!
         *    Calculates the residual vector and Jacobian due to volume heat
         *    source.
         */
        virtual bool volume_heat_source_residual(bool request_jacobian,
                                                 RealVectorX& f,
                                                 RealMatrixX& jac,
                                                 MAST::BoundaryConditionBase& p);

        
        /*!
         *    Calculates the residual vector and Jacobian due to volume heat 
         *    source.
         */
        virtual bool
        volume_heat_source_residual_sensitivity(bool request_jacobian,
                                                RealVectorX& f,
                                                RealMatrixX& jac,
                                                MAST::BoundaryConditionBase& p);
        
        /*!
         *    When \p mass = false, initializes the FEM operator matrix to the
         *    shape functions as
         *    \f[  B = \left[ \begin{array}{c}
         *    {\bf N} \\ {\bf N} \\ {\bf N}
         *    \end{array} \right] \f]
         */
        void _initialize_mass_fem_operator(const unsigned int qp,
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


#endif // __mast__heat_conduction_elem_base__
