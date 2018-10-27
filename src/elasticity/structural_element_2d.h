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

#ifndef __mast__structural_element_2d__
#define __mast__structural_element_2d__

// C++ includes
#include <memory>


// MAST includes
#include "elasticity/bending_structural_element.h"




namespace MAST {
    
    // Forward declerations
    class BendingOperator2D;
    class BoundaryCondition;
    class FEMOperatorMatrix;
    
    
    
    class StructuralElement2D:
    public MAST::BendingStructuralElem {
        
    public:
        StructuralElement2D(MAST::SystemInitialization& sys,
                            MAST::AssemblyBase& assembly,
                            const libMesh::Elem& elem,
                            const MAST::ElementPropertyCardBase& p);
        
        virtual ~StructuralElement2D();
        
        /*!
         *    row dimension of the direct strain matrix, also used for the
         *    bending operator row dimension
         */
        virtual unsigned int n_direct_strain_components() {
            return 3;
        }
        
        /*!
         *    row dimension of the von Karman strain matrix
         */
        virtual unsigned int n_von_karman_strain_components() {
            return 2;
        }
        
        /*!
         *    Calculates the internal residual vector and Jacobian due to
         *    strain energy
         */
        virtual bool internal_residual(bool request_jacobian,
                                       RealVectorX& f,
                                       RealMatrixX& jac);
        
        /*!
         *    Calculates the sensitivity internal residual vector and Jacobian due to
         *    strain energy
         */
        virtual bool internal_residual_sensitivity(const MAST::FunctionBase& p,
                                                   bool request_jacobian,
                                                   RealVectorX& f,
                                                   RealMatrixX& jac);
        /*!
         *   calculates d[J]/d{x} . d{x}/dp
         */
        virtual bool
        internal_residual_jac_dot_state_sensitivity (RealMatrixX& jac);
        
        /*!
         *   calculates the term on side \par s:
         *   \f$ \int_\Gamma a(\delta u, u) v_n ~d\Gamma \f$.
         *
         */
        virtual void
        internal_residual_boundary_velocity(const MAST::FunctionBase& p,
                                            const unsigned int s,
                                            const MAST::FieldFunction<RealVectorX>& vel_f,
                                            bool request_jacobian,
                                            RealVectorX& f,
                                            RealMatrixX& jac);

        /*!
         *    Calculates the internal residual vector and Jacobian due to
         *    strain energy coming from a prestress
         */
        virtual bool prestress_residual (bool request_jacobian,
                                         RealVectorX& f,
                                         RealMatrixX& jac);
        
        /*!
         *    Calculates the internal residual vector and Jacobian due to
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
         *    Calculates the stress tensor
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
                                           const MAST::FieldFunction<RealVectorX>& vel_f);
    

        virtual void
        calculate_stress_temperature_derivative(MAST::FEBase& fe_thermal,
                                                MAST::StressStrainOutputBase& output);
        
        virtual void
        thermal_residual_temperature_derivative (const MAST::FEBase& fe_thermal,
                                                 RealMatrixX& m);
    
    protected:
        
        
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
         *    Calculates the residual vector and Jacobian due to thermal stresses
         */
        virtual bool thermal_residual(bool request_jacobian,
                                      RealVectorX& f,
                                      RealMatrixX& jac,
                                      MAST::BoundaryConditionBase& bc);
        
        /*!
         *    Calculates the sensitivity of residual vector and the Jacobian due to
         *    thermal stresses
         */
        virtual bool thermal_residual_sensitivity(const MAST::FunctionBase& p,
                                                  bool request_jacobian,
                                                  RealVectorX& f,
                                                  RealMatrixX& jac,
                                                  MAST::BoundaryConditionBase& bc);

        /*!
         *   calculates the term on side \par s:
         *   \f$ \int_\Gamma a(\delta u, u) v_n ~d\Gamma \f$.
         *
         */
        virtual void
        thermal_residual_boundary_velocity(const MAST::FunctionBase& p,
                                           const unsigned int s,
                                           const MAST::FieldFunction<RealVectorX>& vel_f,
                                           MAST::BoundaryConditionBase& bc,
                                           bool request_jacobian,
                                           RealVectorX& f,
                                           RealMatrixX& jac);

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
                                            MAST::BoundaryConditionBase& bc);

        
        /*!
         *    Calculates the force vector and Jacobian due to piston-theory
         *    based surface pressure on the side identified by \p side.
         *    The order of the boundary condition and direction of fluid
         *     flow are obtained from the BoundaryConditionBase object.
         */
        virtual bool
        piston_theory_residual_sensitivity(const MAST::FunctionBase& p,
                                           bool request_jacobian,
                                           RealVectorX &f,
                                           RealMatrixX& jac_xdot,
                                           RealMatrixX& jac,
                                           const unsigned int side,
                                           MAST::BoundaryConditionBase& bc);

        
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
         *   initialze the von Karman strain in \par vK_strain, the operator
         *   matrices needed for Jacobian calculation.
         *   vk_strain = [dw/dx 0; 0 dw/dy; dw/dy dw/dx]
         *   Bmat_vk   = [dw/dx; dw/dy]
         */
        virtual void
        initialize_von_karman_strain_operator(const unsigned int qp,
                                              const MAST::FEBase& fe,
                                              RealVectorX& vk_strain,
                                              RealMatrixX& vk_dwdxi_mat,
                                              MAST::FEMOperatorMatrix& Bmat_vk);
        
        /*!
         *   initialze the sensitivity of von Karman operator
         *   matrices needed for Jacobian calculation.
         *   vk_dwdxi_mat_sens = [dw/dx 0; 0 dw/dy; dw/dy dw/dx]
         */
        virtual void
        initialize_von_karman_strain_operator_sensitivity(const unsigned int qp,
                                                          const MAST::FEBase& fe,
                                                          RealMatrixX& vk_dwdxi_mat_sens);
        

        virtual void
        initialize_green_lagrange_strain_operator(const unsigned int qp,
                                                  const MAST::FEBase& fe,
                                                  const RealVectorX& local_disp,
                                                  RealVectorX& epsilon,
                                                  RealMatrixX& mat_x,
                                                  RealMatrixX& mat_y,
                                                  MAST::FEMOperatorMatrix& Bmat_lin,
                                                  MAST::FEMOperatorMatrix& Bmat_nl_x,
                                                  MAST::FEMOperatorMatrix& Bmat_nl_y,
                                                  MAST::FEMOperatorMatrix& Bmat_nl_u,
                                                  MAST::FEMOperatorMatrix& Bmat_nl_v);

        /*!
         *   performs integration at the quadrature point for the provided
         *   matrices. The temperature vector and matrix entities are provided for
         *   integration
         */
        virtual void
        _internal_residual_operation(bool                       if_bending,
                                     bool                       if_vk,
                                     const unsigned int         n2,
                                     const unsigned int         qp,
                                     const MAST::FEBase&        fe,
                                     const std::vector<Real>&   JxW,
                                     bool                       request_jacobian,
                                     RealVectorX&               local_f,
                                     RealMatrixX&               local_jac,
                                     RealVectorX&               local_disp,
                                     RealVectorX&               strain_mem,
                                     FEMOperatorMatrix&         Bmat_lin,
                                     FEMOperatorMatrix&         Bmat_nl_x,
                                     FEMOperatorMatrix&         Bmat_nl_y,
                                     FEMOperatorMatrix&         Bmat_nl_u,
                                     FEMOperatorMatrix&         Bmat_nl_v,
                                     MAST::FEMOperatorMatrix&   Bmat_bend,
                                     MAST::FEMOperatorMatrix&   Bmat_vk,
                                     RealMatrixX&               mat_x,
                                     RealMatrixX&               mat_y,
                                     RealMatrixX&               stress,
                                     RealMatrixX&               vk_dwdxi_mat,
                                     RealMatrixX&               material_A_mat,
                                     RealMatrixX&               material_B_mat,
                                     RealMatrixX&               material_D_mat,
                                     RealVectorX&               vec1_n1,
                                     RealVectorX&               vec2_n1,
                                     RealVectorX&               vec3_n2,
                                     RealVectorX&               vec4_2,
                                     RealVectorX&               vec5_2,
                                     RealVectorX&               vec6_n2,
                                     RealMatrixX&               mat1_n1n2,
                                     RealMatrixX&               mat2_n2n2,
                                     RealMatrixX&               mat3,
                                     RealMatrixX&               mat4_2n2,
                                     RealMatrixX&               mat5_3n2);
        
        
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

        /*!
         *    bending operator used for this elmeent
         */
        MAST::BendingOperator2D *_bending_operator;
    };
}



#endif // __mast__structural_element_2d__
