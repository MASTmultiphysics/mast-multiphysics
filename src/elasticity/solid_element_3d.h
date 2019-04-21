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

#ifndef __mast__solid_element_3d__
#define __mast__solid_element_3d__

/// MAST includes
#include "elasticity/structural_element_base.h"

// Forward declerations
class FEMOperatorMatrix;

namespace MAST {
    
    // Forward declerations
    class BoundaryConditionBase;
    
    
    class StructuralElement3D:
    public MAST::StructuralElementBase {
        
    public:
        StructuralElement3D(MAST::SystemInitialization& sys,
                            MAST::AssemblyBase& assembly,
                            const MAST::GeomElem& elem,
                            const MAST::ElementPropertyCardBase& p);
        
        
        virtual ~StructuralElement3D();
        
        /*!
         *   Calculates the inertial force and the Jacobian matrices
         */
        virtual bool inertial_residual (bool request_jacobian,
                                        RealVectorX& f,
                                        RealMatrixX& jac_xddot,
                                        RealMatrixX& jac_xdot,
                                        RealMatrixX& jac);
        
        
        /*!
         *    Calculates the internal residual vector and Jacobian due to
         *    strain energy
         */
        virtual bool internal_residual(bool request_jacobian,
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
         *    Calculates the sensitivity of internal residual vector and
         *    Jacobian due to strain energy
         */
        virtual bool internal_residual_sensitivity(const MAST::FunctionBase& p,
                                                   bool request_jacobian,
                                                   RealVectorX& f,
                                                   RealMatrixX& jac);
        
        /*!
         *   calculates d[J]/d{x} . d{x}/dp
         */
        virtual bool
        internal_residual_jac_dot_state_sensitivity (RealMatrixX& jac) {
            libmesh_assert(false); // to be implemented for 3D elements
        }
        
        /*!
         *    Calculates the prestress residual vector and Jacobian
         */
        virtual bool prestress_residual (bool request_jacobian,
                                         RealVectorX& f,
                                         RealMatrixX& jac);
        
        
        /*!
         *    Calculates the sensitivity prestress residual vector and Jacobian
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
         *  @returns the dimension of the incompatible mode vector
         */
        virtual unsigned int incompatible_mode_size() const {
            return 30;
        }

        /*!
         *    updates the incompatible solution for this element. \p dsol
         *    is the update to the element solution for the current
         *    nonlinear step.
         */
        virtual void update_incompatible_mode_solution(const RealVectorX& dsol);
         
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
         *    Calculates the sensitivity of residual vector and Jacobian due to
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
                                            MAST::BoundaryConditionBase& bc) {
            libmesh_error_msg("Invalid Function Call: Piston theory force \
                              is not defined for solid element volume.");
        }

        
        /*!
         *    Calculates the force vector and Jacobian sensitivity due to 
         *    piston-theory based surface pressure on the entire element domain.
         *    This is applicable for only 1D and 2D elements. The order
         *    of the boundary condition and direction of fluid flow are
         *    obtained from the BoundaryConditionBase object.
         */
        virtual bool piston_theory_residual_sensitivity(const MAST::FunctionBase& p,
                                                        bool request_jacobian,
                                                        RealVectorX &f,
                                                        RealMatrixX& jac_xdot,
                                                        RealMatrixX& jac,
                                                        MAST::BoundaryConditionBase& bc) {
            
            libmesh_error_msg("Invalid Function Call: Piston theory force \
                              is not defined for solid element volume.");
        }

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
         *    Calculates the force vector and Jacobian sensitivity  due to
         *    piston-theory based surface pressure on the side identified by
         *     \p side. The order of the boundary condition and direction of
         *     fluid flow are obtained from the BoundaryConditionBase object.
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
            return false;
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
                                           const MAST::FieldFunction<RealVectorX>& vel_f) {
            
            libmesh_error(); // to be implemented
        }

        /*!
         *   initialize strain operator matrix
         */
        void initialize_strain_operator (const unsigned int qp,
                                         const MAST::FEBase& fe,
                                         FEMOperatorMatrix& Bmat);
        
        /*!
         *    initialize the strain operator matrices for the 
         *    Green-Lagrange strain matrices
         */
        void initialize_green_lagrange_strain_operator(const unsigned int qp,
                                                       const MAST::FEBase& fe,
                                                       const RealVectorX& local_disp,
                                                       RealVectorX& epsilon,
                                                       RealMatrixX& mat_x,
                                                       RealMatrixX& mat_y,
                                                       RealMatrixX& mat_z,
                                                       MAST::FEMOperatorMatrix& Bmat_lin,
                                                       MAST::FEMOperatorMatrix& Bmat_nl_x,
                                                       MAST::FEMOperatorMatrix& Bmat_nl_y,
                                                       MAST::FEMOperatorMatrix& Bmat_nl_z,
                                                       MAST::FEMOperatorMatrix& Bmat_nl_u,
                                                       MAST::FEMOperatorMatrix& Bmat_nl_v,
                                                       MAST::FEMOperatorMatrix& Bmat_nl_w);
        
        /*!
         *   initialize incompatible strain operator
         */
        void initialize_incompatible_strain_operator(const unsigned int qp,
                                                     const MAST::FEBase& fe,
                                                     FEMOperatorMatrix& Bmat,
                                                     RealMatrixX& G_mat);

        /*!
         *   initialize the Jacobian needed for incompatible modes
         */
        void _init_incompatible_fe_mapping( const libMesh::Elem& e);

        
        /*!
         *   Jacobian matrix at element center needed for incompatible modes
         */
        RealMatrixX _T0_inv_tr;

    };
}

#endif // __mast__solid_element_3d__
