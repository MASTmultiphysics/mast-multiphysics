
#ifndef __mast__structural_element_2d__
#define __mast__structural_element_2d__

// C++ includes
#include <memory>


// MAST includes
#include "elasticity/bending_structural_element.h"




namespace MAST {
    
    // Forward declerations
    class BoundaryCondition;
    class FEMOperatorMatrix;
    
    
    
    class StructuralElement2D:
    public MAST::BendingStructuralElem {
        
    public:
        StructuralElement2D(MAST::SystemInitialization& sys,
                            const libMesh::Elem& elem,
                            const MAST::ElementPropertyCardBase& p);
        
        
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
                                       RealMatrixX& jac,
                                       bool if_ignore_ho_jac);
        
        /*!
         *    Calculates the sensitivity internal residual vector and Jacobian due to
         *    strain energy
         */
        virtual bool internal_residual_sensitivity(bool request_jacobian,
                                                   RealVectorX& f,
                                                   RealMatrixX& jac,
                                                   bool if_ignore_ho_jac);
        /*!
         *   calculates d[J]/d{x} . d{x}/dp
         */
        virtual bool
        internal_residual_jac_dot_state_sensitivity (RealMatrixX& jac) {
            libmesh_assert(false); // to be implemented for 2D elements
        }
        
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
        virtual bool prestress_residual_sensitivity (bool request_jacobian,
                                                     RealVectorX& f,
                                                     RealMatrixX& jac);
        
    protected:
        
        
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
        virtual bool thermal_residual_sensitivity(bool request_jacobian,
                                                  RealVectorX& f,
                                                  RealMatrixX& jac,
                                                  MAST::BoundaryConditionBase& bc);
        
        /*!
         *   initialize membrane strain operator matrix
         */
        virtual void
        initialize_direct_strain_operator(const unsigned int qp,
                                          MAST::FEMOperatorMatrix& Bmat);
        
        /*!
         *   initialze the von Karman strain in \par vK_strain, the operator
         *   matrices needed for Jacobian calculation.
         *   vk_strain = [dw/dx 0; 0 dw/dy; dw/dy dw/dx]
         *   Bmat_vk   = [dw/dx; dw/dy]
         */
        virtual void
        initialize_von_karman_strain_operator(const unsigned int qp,
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
                                                          RealMatrixX& vk_dwdxi_mat_sens);
        
        /*!
         *   performs integration at the quadrature point for the provided
         *   matrices. The temperature vector and matrix entities are provided for
         *   integration
         */
        virtual void
        _internal_residual_operation(bool if_bending,
                                     bool if_vk,
                                     const unsigned int n2,
                                     const unsigned int qp,
                                     const std::vector<Real>& JxW,
                                     bool request_jacobian,
                                     bool if_ignore_ho_jac,
                                     RealVectorX& local_f,
                                     RealMatrixX& local_jac,
                                     MAST::FEMOperatorMatrix& Bmat_mem,
                                     MAST::FEMOperatorMatrix& Bmat_bend,
                                     MAST::FEMOperatorMatrix& Bmat_vk,
                                     RealMatrixX& stress,
                                     RealMatrixX& stress_l,
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
         *   sensitivity of linear part of the geometric stiffness matrix
         */
        void
        _linearized_geometric_stiffness_sensitivity_with_static_solution
        (const unsigned int n2,
         const unsigned int qp,
         const std::vector<Real>& JxW,
         RealMatrixX& local_jac,
         MAST::FEMOperatorMatrix& Bmat_mem,
         MAST::FEMOperatorMatrix& Bmat_bend,
         MAST::FEMOperatorMatrix& Bmat_vk,
         RealMatrixX& stress_l,
         RealMatrixX& vk_dwdxi_mat,
         RealMatrixX& material_A_mat,
         RealMatrixX& material_B_mat,
         RealVectorX& vec1_n1,
         RealVectorX& vec2_n1,
         RealMatrixX& mat1_n1n2,
         RealMatrixX& mat2_n2n2,
         RealMatrixX& mat3);
        
        
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



#endif // __mast__structural_element_2d__
