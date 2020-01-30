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

#ifndef __mast__structural_element_1d__
#define __mast__structural_element_1d__

// MAST includes
#include "elasticity/bending_structural_element.h"




namespace MAST {
    
    // Forward declerations
    class BendingOperator1D;
    class BoundaryConditionBase;
    class FEMOperatorMatrix;
    
    
    class StructuralElement1D:
    public MAST::BendingStructuralElem {
        
    public:
        StructuralElement1D(MAST::SystemInitialization& sys,
                            const MAST::GeomElem& elem,
                            const MAST::ElementPropertyCardBase& p);
        
        
        virtual ~StructuralElement1D();
        
        
        /*!
         *    row dimension of the direct strain matrix, also used for the
         *    bending operator row dimension
         */
        virtual unsigned int n_direct_strain_components() {
            return 2;
        }
        
        /*!
         *    row dimension of the von Karman strain matrix
         */
        virtual unsigned int n_von_karman_strain_components() {
            return 2;
        }
        
        /*!
         *    Calculates the internal force vector and Jacobian due to
         *    strain energy
         */
        virtual bool internal_residual(bool request_jacobian,
                                       RealVectorX& f,
                                       RealMatrixX& jac);
        
        /*!
         *    Calculates the sensitivity internal force vector and Jacobian due to
         *    strain energy
         */
        virtual bool internal_residual_sensitivity(const MAST::FunctionBase& p,
                                                   bool request_jacobian,
                                                   RealVectorX& f,
                                                   RealMatrixX& jac);
        
        /*!
         *   calculates the term on side \p s:
         *   \f$ \int_\Gamma a(\delta u, u) v_n ~d\Gamma \f$.
         *
         */
        virtual void
        internal_residual_boundary_velocity(const MAST::FunctionBase& p,
                                            const unsigned int s,
                                            const MAST::FieldFunction<RealVectorX>& vel_f,
                                            bool request_jacobian,
                                            RealVectorX& f,
                                            RealMatrixX& jac) {
            libmesh_assert(false);
        }

        /*!
         *   calculates d[J]/d{x} . d{x}/dp
         */
        virtual bool
        internal_residual_jac_dot_state_sensitivity (RealMatrixX& jac);
        
        /*!
         *    Calculates the internal force vector and Jacobian due to
         *    strain energy coming from a prestress
         */
        virtual bool prestress_residual (bool request_jacobian,
                                         RealVectorX& f,
                                         RealMatrixX& jac);
        
        /*!
         *    Calculates the internal force vector and Jacobian due to
         *    strain energy coming from a prestress
         */
        virtual bool prestress_residual_sensitivity (const MAST::FunctionBase& p,
                                                     bool request_jacobian,
                                                     RealVectorX& f,
                                                     RealMatrixX& jac);
        
        /*!
         *  @returns false since this element formulation does not use
         *  incompatible modes
         */
        virtual bool if_incompatible_modes() const {
            return false;
        }

        /*!
         *    Calculates the stress tensor. If derivative and sensitivity
         *    with respect to the parameter \p sesitivity_param are calculated
         *    and provided if the respective flags are true.
         */
        virtual bool calculate_stress(bool request_derivative,
                                      const MAST::FunctionBase* p,
                                      MAST::StressStrainOutputBase& output);
        
        /*!
         *    Calculates the boundary velocity term contributions to the
         *    sensitivity of stress at the specified boundary of this element.
         */
        virtual void
        calculate_stress_boundary_velocity(const MAST::FunctionBase& p,
                                           MAST::StressStrainOutputBase& output,
                                           const unsigned int s,
                                           const MAST::FieldFunction<RealVectorX>& vel_f) {
            libmesh_error(); // to be implemented
        }

        virtual void
        calculate_stress_temperature_derivative(MAST::FEBase& fe_thermal,
                                                MAST::StressStrainOutputBase& output) {
            libmesh_error(); // to be implemented
        }
        
        virtual void
        thermal_residual_temperature_derivative (const MAST::FEBase& fe_thermal,
                                                 RealMatrixX& m) {
            libmesh_error(); // to be implemented
        }

    protected:
        
        /*!
         *    Calculates the force vector and Jacobian due to surface traction.
         */
        virtual bool
        surface_traction_residual(bool request_jacobian,
                                  RealVectorX& f,
                                  RealMatrixX& jac,
                                  const unsigned int side,
                                  MAST::BoundaryConditionBase& bc){ return false;}


        /*!
         *    Calculates the sensitivity of element vector and matrix quantities for surface traction
         *    boundary condition.
         */
        virtual bool
        surface_traction_residual_sensitivity(const MAST::FunctionBase& p,
                                              bool request_jacobian,
                                              RealVectorX& f,
                                              RealMatrixX& jac,
                                              const unsigned int side,
                                              MAST::BoundaryConditionBase& bc) {}


        /*!
         *    Calculates the sensitivity of force vector and Jacobian due to surface traction and
         *    sensitiity due to boundary movement.
         */
        virtual bool
        surface_traction_residual_shifted_boundary(bool request_jacobian,
                                                   RealVectorX& f,
                                                   RealMatrixX& jac,
                                                   const unsigned int side,
                                                   MAST::BoundaryConditionBase& bc) {}

        /*!
         *    Calculates the sensitivity of force vector and Jacobian due to surface traction and
         *    sensitiity due to boundary movement.
         */
        virtual bool
        surface_traction_residual_shifted_boundary_sensitivity(const MAST::FunctionBase& p,
                                                               bool request_jacobian,
                                                               RealVectorX& f,
                                                               RealMatrixX& jac,
                                                               const unsigned int side,
                                                               MAST::BoundaryConditionBase& bc) {}



        /*!
         *    Calculates the force vector and Jacobian due to surface pressure.
         */
        bool surface_pressure_residual(bool request_jacobian,
                                       RealVectorX &f,
                                       RealMatrixX &jac,
                                       const unsigned int side,
                                       MAST::BoundaryConditionBase& bc);
        
        /*!
         *    Calculates the force vector and Jacobian due to surface pressure.
         */
        virtual bool
        surface_pressure_residual_sensitivity(const MAST::FunctionBase& p,
                                              bool request_jacobian,
                                              RealVectorX& f,
                                              RealMatrixX& jac,
                                              const unsigned int side,
                                              MAST::BoundaryConditionBase& bc);

        
        /*!
         *    Calculates the force vector and Jacobian due to thermal stresses
         */
        virtual bool thermal_residual(bool request_jacobian,
                                      RealVectorX& f,
                                      RealMatrixX& jac,
                                      MAST::BoundaryConditionBase& bc);
        
        /*!
         *    Calculates the sensitivity of force vector and the Jacobian due to
         *    thermal stresses
         */
        virtual bool thermal_residual_sensitivity(const MAST::FunctionBase& p,
                                                  bool request_jacobian,
                                                  RealVectorX& f,
                                                  RealMatrixX& jac,
                                                  MAST::BoundaryConditionBase& bc);
        
        
        /*!
         *    Calculates the sensitivity of force vector and Jacobian due to
         *    thermal stresses. this should be implemented for each element type
         */
        virtual void thermal_residual_boundary_velocity(const MAST::FunctionBase& p,
                                                        const unsigned int s,
                                                        const MAST::FieldFunction<RealVectorX>& vel_f,
                                                        MAST::BoundaryConditionBase& bc,
                                                        bool request_jacobian,
                                                        RealVectorX& f,
                                                        RealMatrixX& jac) {
            
            libmesh_assert(false); // to be implemented
        }

        /*!
         *    Calculates the force vector and Jacobian due to piston-theory
         *    based surface pressure on the entire element domain.
         *    This is applicable for only 1D and 2D elements. The order
         *    of the boundary condition and direction of fluid flow are
         *    obtained from the BoundaryConditionBase object.
         */
        virtual bool piston_theory_residual(bool request_jacobian,
                                            RealVectorX &f,
                                            RealMatrixX& jac_xdot,
                                            RealMatrixX& jac,
                                            MAST::BoundaryConditionBase& bc);

        
        /*!
         *    Calculates the sensitivity of force vector and
         *    Jacobian due to piston-theory based surface pressure on the
         *    entire element domain. This is applicable for only 1D and 2D
         *    elements. The order of the boundary condition and direction of
         *    fluid flow are obtained from the BoundaryConditionBase object.
         */
        virtual bool
        piston_theory_residual_sensitivity(const MAST::FunctionBase& p,
                                           bool request_jacobian,
                                           RealVectorX &f,
                                           RealMatrixX& jac_xdot,
                                           RealMatrixX& jac,
                                           MAST::BoundaryConditionBase& bc);

        
        /*!
         *    Calculates the force vector and Jacobian due to piston-theory
         *    based surface pressure on the side identified by \p side.
         *    The order of the boundary condition and direction of fluid
         *     flow are obtained from the BoundaryConditionBase object.
         */
        virtual bool piston_theory_residual(bool request_jacobian,
                                            RealVectorX &f,
                                            RealMatrixX& jac_xdot,
                                            RealMatrixX& jac,
                                            const unsigned int side,
                                            MAST::BoundaryConditionBase& bc) {
            
            libmesh_error_msg("Invalid Function Call: Piston theory force \
                              is not defined for 1D element side.");
        }

        
        /*!
         *    Calculates the force vector and Jacobian due to piston-theory
         *    sensitivity based surface pressure on the side identified by \p side.
         *    The order of the boundary condition and direction of fluid
         *     flow are obtained from the BoundaryConditionBase object.
         */
        virtual bool piston_theory_residual_sensitivity(const MAST::FunctionBase& p,
                                                        bool request_jacobian,
                                                        RealVectorX &f,
                                                        RealMatrixX& jac_xdot,
                                                        RealMatrixX& jac,
                                                        const unsigned int side,
                                                        MAST::BoundaryConditionBase& bc) {
            
            libmesh_error_msg("Invalid Function Call: Piston theory force sensitivity \
                              is not defined for 1D element side.");
        }


        /*!
         *    Calculates the force vector and Jacobian due to small
         *    perturbation surface pressure.
         */
        virtual bool
        linearized_frequency_domain_surface_pressure_residual
        (bool request_jacobian,
         ComplexVectorX& f,
         ComplexMatrixX& jac,
         const unsigned int side,
         MAST::BoundaryConditionBase& bc) {
            
            libmesh_error(); // to be implemented
        }
        
        
        
        /*!
         *    Calculates the sensitivity of force vector and Jacobian due to small
         *     is applicable for perturbation surface pressure.
         */
        virtual bool
        linearized_frequency_domain_surface_pressure_residual_sensitivity
        (const MAST::FunctionBase& p,
         bool request_jacobian,
         ComplexVectorX& f,
         ComplexMatrixX& jac,
         const unsigned int side,
         MAST::BoundaryConditionBase& bc) {
            
            libmesh_error(); // to be implemented
        }
        

        /*!
         *   initialize membrane strain operator matrix
         */
        virtual void initialize_direct_strain_operator(const unsigned int qp,
                                                       const MAST::FEBase& fe,
                                                       MAST::FEMOperatorMatrix& Bmat);
        
        /*!
         *   initialze the von Karman strain in \p vK_strain, the operator
         *   matrices needed for Jacobian calculation.
         *   vk_strain = [dw/dx 0; 0 dw/dy; dw/dy dw/dx]
         *   Bmat_vk   = [dw/dx; dw/dy]
         */
        virtual void
        initialize_von_karman_strain_operator(const unsigned int qp,
                                              const MAST::FEBase& fe,
                                              RealVectorX& vk_strain,
                                              RealMatrixX& vk_dvdxi_mat,
                                              RealMatrixX& vk_dwdxi_mat,
                                              MAST::FEMOperatorMatrix& Bmat_v_vk,
                                              MAST::FEMOperatorMatrix& Bmat_w_vk);
        
        /*!
         *   initialze the sensitivity of von Karman operator
         *   matrices needed for Jacobian calculation.
         *   vk_dwdxi_mat_sens = [dw/dx 0; 0 dw/dy; dw/dy dw/dx]
         */
        void
        initialize_von_karman_strain_operator_sensitivity
        (const unsigned int qp,
         const MAST::FEBase& fe,
         RealMatrixX& vk_dvdxi_mat_sens,
         RealMatrixX& vk_dwdxi_mat_sens);
        
        
        /*!
         *   performs integration at the quadrature point for the provided
         *   matrices. The temperature vector and matrix entities are provided for
         *   integration
         */
        virtual void _internal_residual_operation(bool if_vk,
                                                  const unsigned int n2,
                                                  const unsigned int qp,
                                                  const MAST::FEBase& fe,
                                                  const std::vector<Real>& JxW,
                                                  bool request_jacobian,
                                                  RealVectorX& local_f,
                                                  RealMatrixX& local_jac,
                                                  MAST::BendingOperator1D* bend_op,
                                                  MAST::FEMOperatorMatrix& Bmat_mem,
                                                  MAST::FEMOperatorMatrix& Bmat_bend_v,
                                                  MAST::FEMOperatorMatrix& Bmat_bend_w,
                                                  MAST::FEMOperatorMatrix& Bmat_v_vk,
                                                  MAST::FEMOperatorMatrix& Bmat_w_vk,
                                                  RealMatrixX& stress,
                                                  RealMatrixX& stress_l,
                                                  RealMatrixX& vk_dvdxi_mat,
                                                  RealMatrixX& vk_dwdxi_mat,
                                                  RealMatrixX& material_A_mat,
                                                  RealMatrixX& material_B_mat,
                                                  RealMatrixX& material_D_mat,
                                                  RealVectorX& vec1_n1,
                                                  RealVectorX& vec2_n1,
                                                  RealVectorX& vec3_n2,
                                                  RealVectorX& vec4_2,
                                                  RealVectorX& vec5_2,
                                                  RealMatrixX& mat1_n1n2,
                                                  RealMatrixX& mat2_n2n2,
                                                  RealMatrixX& mat3,
                                                  RealMatrixX& mat4_2n2);
        
        

        /*!
         *   converts the prestress stress tensor to a vector representation
         */
        void _convert_prestress_A_mat_to_vector(const RealMatrixX& mat,
                                                RealVectorX& vec) const;

        /*!
         *   converts the prestress stress tensor to a vector representation
         */
        void _convert_prestress_B_mat_to_vector(const RealMatrixX& mat,
                                                RealVectorX& vec) const;

    };
}


#endif // __mast__structural_element_1d__
