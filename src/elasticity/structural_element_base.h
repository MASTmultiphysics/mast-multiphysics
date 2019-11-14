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

#ifndef __mast__structural_element_base__
#define __mast__structural_element_base__

// C++ includes
#include <memory>
#include <map>

// MAST includes
#include "base/elem_base.h"


namespace MAST {
    
    // Forward declerations
    class GeomElem;
    class ElementPropertyCardBase;
    class BoundaryConditionBase;
    class FEMOperatorMatrix;
    class StressStrainOutputBase;
    template <typename ValType> class FieldFunction;
    
    
    class StructuralElementBase:
    public MAST::ElementBase
    {
    public:
        /*!
         *   Constructor.
         */
        StructuralElementBase(MAST::SystemInitialization& sys,
                              MAST::AssemblyBase& assembly,
                              const MAST::GeomElem& elem,
                              const MAST::ElementPropertyCardBase& p);
        
        virtual ~StructuralElementBase();

        
        /*!
         *   stores \p vec as solution for element level calculations,
         *   or its sensitivity if \p if_sens is true.
         */
        virtual void set_solution(const RealVectorX& vec,
                                  bool if_sens = false);

        
        /*!
         *   stores \p vec as perturbed solution for element level
         *   calculations, or its sensitivity if \p if_sens is true.
         */
        virtual void set_perturbed_solution(const RealVectorX& vec,
                                            bool if_sens = false);

        
        /*!
         *    stores \p vec as velocity for element level calculations,
         *    or its sensitivity if \p if_sens is true.
         */
        virtual void set_velocity(const RealVectorX& vec,
                                  bool if_sens = false);

        /*!
         *    stores \p vec as perturbed velocity for element level
         *    calculations, or its sensitivity if \p if_sens is true.
         */
        virtual void set_perturbed_velocity(const RealVectorX& vec,
                                            bool if_sens = false);

        
        /*!
         *    stores \p vec as acceleration for element level
         *    calculations, or its sensitivity if \p if_sens is true.
         */
        virtual void set_acceleration(const RealVectorX& vec,
                                      bool if_sens = false);

        
        /*!
         *    stores \p vec as perturbed acceleration for element level
         *    calculations, or its sensitivity if \p if_sens is true.
         */
        virtual void set_perturbed_acceleration(const RealVectorX& vec,
                                                bool if_sens = false);

        
        /*!
         *  @returns a constant reference to the element solution 
         *  (or its derivative if \p if_sens is true) in the local
         *  element coordinate system
         */
        const RealVectorX& local_solution(bool if_sens = false) const {
            
            if (!if_sens)
                return _local_sol;
            else
                return _local_sol_sens;
        }
        
        /*!
         *   returns a constant reference to the finite element object
         */
        const MAST::ElementPropertyCardBase& elem_property()  {
            return _property;
        }
        
        
        /*!
         *   internal force contribution to system residual
         */
        virtual bool internal_residual (bool request_jacobian,
                                        RealVectorX& f,
                                        RealMatrixX& jac) = 0;

        /*!
         *   internal force contribution to system residual of the linearized
         *   problem
         */
        virtual bool
        linearized_internal_residual (bool request_jacobian,
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
                                            RealMatrixX& jac) = 0;
        
        /*!
         *   calculates d[J]/d{x} . d{x}/dp
         */
        virtual bool
        internal_residual_jac_dot_state_sensitivity (RealMatrixX& jac) = 0;
        
        
        
        /*!
         *   damping force contribution to system residual
         */
        virtual bool damping_residual (bool request_jacobian,
                                       RealVectorX& f,
                                       RealMatrixX& jac_xdot,
                                       RealMatrixX& jac) {
            
            // to be implemented nothing is done yet
            return false;
        }
        
        /*!
         *   inertial force contribution to system residual
         */
        virtual bool inertial_residual (bool request_jacobian,
                                        RealVectorX& f,
                                        RealMatrixX& jac_xddot,
                                        RealMatrixX& jac_xdot,
                                        RealMatrixX& jac);

        /*!
         *   inertial force contribution to system residual of linerized 
         *   problem
         */
        virtual bool linearized_inertial_residual (bool request_jacobian,
                                                   RealVectorX& f,
                                                   RealMatrixX& jac_xddot,
                                                   RealMatrixX& jac_xdot,
                                                   RealMatrixX& jac);

        
        /*!
         *   side external force contribution to system residual. This primarily
         *   handles the boundary conditions. If requested, the Jacobians of the
         *   residual due to xdot will be returned in \p jac_xdot and the 
         *   Jacobian due to x is returned in jac.
         */
        bool side_external_residual (bool request_jacobian,
                                     RealVectorX& f,
                                     RealMatrixX& jac_xdot,
                                     RealMatrixX& jac,
                                     std::multimap<libMesh::boundary_id_type, MAST::BoundaryConditionBase*>& bc);

        /*!
         *   side external force contribution to system residual. This primarily
         *   handles the boundary conditions. If requested, the Jacobians of the
         *   residual due to xdot will be returned in \p jac_xdot and the
         *   Jacobian due to x is returned in jac.
         */
        bool
        linearized_side_external_residual
        (bool request_jacobian,
         RealVectorX& f,
         RealMatrixX& jac_xdot,
         RealMatrixX& jac,
         std::multimap<libMesh::boundary_id_type, MAST::BoundaryConditionBase*>& bc);
        
        

        /*!
         *   Calculates the external force due to frequency domain
         *   side external force contribution to system residual. This primarily
         *   handles the boundary conditions. If requested, the
         *   Jacobian due to x is returned in jac.
         */
        bool
        linearized_frequency_domain_side_external_residual
        (bool request_jacobian,
         ComplexVectorX& f,
         ComplexMatrixX& jac,
         std::multimap<libMesh::boundary_id_type, MAST::BoundaryConditionBase*>& bc);
        
        
        /*!
         *   prestress force contribution to system residual
         */
        virtual bool prestress_residual (bool request_jacobian,
                                         RealVectorX& f,
                                         RealMatrixX& jac) = 0;
        
        /*!
         *   volume external force contribution to system residual. If 
         *   requested, the Jacobians of the residual due to xdot will be 
         *   returned in \p jac_xdot and the Jacobian due to x is
         *   returned in jac.
         */
        bool volume_external_residual (bool request_jacobian,
                                       RealVectorX& f,
                                       RealMatrixX& jac_xdot,
                                       RealMatrixX& jac,
                                       std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*>& bc);

        
        /*!
         *   boundary velocity contribution of volume external force.
         */
        void volume_external_residual_boundary_velocity
        (const MAST::FunctionBase& p,
         const unsigned int s,
         const MAST::FieldFunction<RealVectorX>& vel_f,
         std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*>& bc,
         bool request_jacobian,
         RealVectorX& f,
         RealMatrixX& jac);

        
        /*!
         *   volume external force contribution to system residual. If
         *   requested, the Jacobians of the residual due to xdot will be
         *   returned in \p jac_xdot and the Jacobian due to x is
         *   returned in jac.
         */
        bool
        linearized_volume_external_residual
        (bool request_jacobian,
         RealVectorX& f,
         RealMatrixX& jac_xdot,
         RealMatrixX& jac,
         std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*>& bc);

        
        /*!
         *   Calculates the frequency domain  volume external force contribution
         *   to system residual. If requested, the Jacobian due to x is
         *   returned in jac.
         */
        bool
        linearized_frequency_domain_volume_external_residual
        (bool request_jacobian,
         ComplexVectorX& f,
         ComplexMatrixX& jac,
         std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*>& bc);
        
        
        
        /*!
         *   sensitivity of the internal force contribution to system residual
         */
        virtual bool internal_residual_sensitivity (const MAST::FunctionBase& p,
                                                    bool request_jacobian,
                                                    RealVectorX& f,
                                                    RealMatrixX& jac) = 0;

        /*!
         *   sensitivity of the damping force contribution to system residual
         */
        virtual bool damping_residual_sensitivity (const MAST::FunctionBase& p,
                                                   bool request_jacobian,
                                                   RealVectorX& f,
                                                   RealMatrixX& jac_xdot,
                                                   RealMatrixX& jac) {
            return false;
        }
        
        /*!
         *   sensitivity of the inertial force contribution to system residual
         */
        virtual bool inertial_residual_sensitivity (const MAST::FunctionBase& p,
                                                    bool request_jacobian,
                                                    RealVectorX& f,
                                                    RealMatrixX& jac_xddot,
                                                    RealMatrixX& jac_xdot,
                                                    RealMatrixX& jac);

        /*!
         *   sensitivity of the inertial force contribution to system residual
         */
        virtual void
        inertial_residual_boundary_velocity (const MAST::FunctionBase& p,
                                             const unsigned int s,
                                             const MAST::FieldFunction<RealVectorX>& vel_f,
                                             bool request_jacobian,
                                             RealVectorX& f,
                                             RealMatrixX& jac_xddot,
                                             RealMatrixX& jac_xdot,
                                             RealMatrixX& jac);

        /*!
         *   sensitivity of the side external force contribution to system residual
         */
        bool side_external_residual_sensitivity (const MAST::FunctionBase& p,
                                                 bool request_jacobian,
                                                 RealVectorX& f,
                                                 RealMatrixX& jac_xdot,
                                                 RealMatrixX& jac,
                                                 std::multimap<libMesh::boundary_id_type, MAST::BoundaryConditionBase*>& bc);
        
        /*!
         *   sensitivity of the prestress force contribution to system residual
         */
        virtual bool prestress_residual_sensitivity (const MAST::FunctionBase& p,
                                                     bool request_jacobian,
                                                     RealVectorX& f,
                                                     RealMatrixX& jac) = 0;
        
        /*!
         *   sensitivity of the volume external force contribution to system residual
         */
        bool volume_external_residual_sensitivity (const MAST::FunctionBase& p,
                                                   bool request_jacobian,
                                                   RealVectorX& f,
                                                   RealMatrixX& jac_xdot,
                                                   RealMatrixX& jac,
                                                   std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*>& bc);
        


        /*!
         *  @returns true/false based on whether or not an element uses
         *  the incompatible mode approach to enrich local strain fields
         */
        virtual bool if_incompatible_modes() const = 0;

        
        /*!
         *  @returns the dimension of the incompatible mode vector
         */
        virtual unsigned int incompatible_mode_size() const {
            // needs to be implemented only in inherited elements
            // that use incompatible modes
            libmesh_error();
            return 0;
        }

        
        /*!
         *  sets the pointer to the incompatible mode solution vector. This
         *  is stored as a pointer, and any modifications will overwrite the
         *  original vector
         */
        void set_incompatible_mode_solution(RealVectorX& vec) {
            _incompatible_sol = &vec;
        }
        
        
        /*!
         *    updates the incompatible solution for this element. \p dsol
         *    is the update to the element solution for the current
         *    nonlinear step.
         */
        virtual void update_incompatible_mode_solution(const RealVectorX& dsol) {
            // should be implemented in derived classes, if relevant
            libmesh_error();
        }

        
        /*!
         *    Calculates the stress tensor. If derivative and sensitivity
         *    with respect to the parameter \p sesitivity_param are calculated
         *    and provided if the respective flags are true.
         */
        virtual bool calculate_stress(bool request_derivative,
                                      const MAST::FunctionBase* f,
                                      MAST::StressStrainOutputBase& output) = 0;
        

        /*!
         *    Calculates the boundary velocity term contributions to the
         *    sensitivity of stress at the specified boundary of this element.
         */
        virtual void
        calculate_stress_boundary_velocity(const MAST::FunctionBase& p,
                                           MAST::StressStrainOutputBase& output,
                                           const unsigned int s,
                                           const MAST::FieldFunction<RealVectorX>& vel_f) = 0;

        virtual void
        calculate_stress_temperature_derivative(MAST::FEBase& fe_thermal,
                                                MAST::StressStrainOutputBase& output) = 0;
        
        virtual void
        thermal_residual_temperature_derivative (const MAST::FEBase& fe_thermal,
                                                 RealMatrixX& m) = 0;

        /*!
         *    flag for follower forces
         */
        bool follower_forces;
        

        template <typename ValType>
        void transform_matrix_to_global_system(const ValType& local_mat,
                                               ValType& global_mat) const;
        
        template <typename ValType>
        void transform_vector_to_local_system(const ValType& global_vec,
                                              ValType& local_vec) const;
        
        template <typename ValType>
        void transform_vector_to_global_system(const ValType& local_vec,
                                               ValType& global_vec) const;
        
    protected:
        
        /*!
         *    Calculates the force vector and Jacobian due to surface traction.
         */
        virtual bool
        surface_traction_residual(bool request_jacobian,
                                  RealVectorX& f,
                                  RealMatrixX& jac,
                                  const unsigned int side,
                                  MAST::BoundaryConditionBase& bc) = 0;

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
                                              MAST::BoundaryConditionBase& bc) = 0;


        /*!
         *    Calculates the sensitivity of force vector and Jacobian due to surface traction and
         *    sensitiity due to boundary movement.
         */
        virtual bool
        surface_traction_residual_shifted_boundary(bool request_jacobian,
                                                   RealVectorX& f,
                                                   RealMatrixX& jac,
                                                   const unsigned int side,
                                                   MAST::BoundaryConditionBase& bc) = 0;

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
                                                               MAST::BoundaryConditionBase& bc) = 0;

        /*!
         *    Calculates the force vector and Jacobian due to surface pressure.
         */
        virtual bool
        surface_pressure_residual(bool request_jacobian,
                                  RealVectorX& f,
                                  RealMatrixX& jac,
                                  const unsigned int side,
                                  MAST::BoundaryConditionBase& bc) = 0;
        
        
        /*!
         *    Calculates the force vector and Jacobian due to surface pressure
         *    applied on the entire element domain. This is applicable for
         *    only 1D and 2D elements.
         */
        virtual bool
        surface_pressure_residual(bool request_jacobian,
                                  RealVectorX& f,
                                  RealMatrixX& jac,
                                  MAST::BoundaryConditionBase& bc);

        /*!
         *    Calculates the force vector and Jacobian due to surface pressure
         *    applied on the entire element domain. This is applicable for
         *    only 1D and 2D elements.
         */
        virtual void
        surface_pressure_boundary_velocity(const MAST::FunctionBase& p,
                                           const unsigned int s,
                                           const MAST::FieldFunction<RealVectorX>& vel_f,
                                           MAST::BoundaryConditionBase& bc,
                                           bool request_jacobian,
                                           RealVectorX& f,
                                           RealMatrixX& jac);

        /*!
         *    Calculates the force vector and Jacobian due to surface pressure
         *    applied on the entire element domain. This is applicable for
         *    only 1D and 2D elements.
         */
        virtual bool
        linearized_surface_pressure_residual(bool request_jacobian,
                                             RealVectorX& f,
                                             RealMatrixX& jac,
                                             MAST::BoundaryConditionBase& bc);

        
        /*!
         *   @returns the piston theory cp value based on the 
         *   specified normal velocity and the flow parameters
         */
        Real piston_theory_cp(const unsigned int order,
                              const Real vel_U,
                              const Real gamma,
                              const Real mach);

        
        
        /*!
         *   @returns the derivative of piston theory cp value with respect 
         *   to the normal velocity. 
         */
        Real piston_theory_dcp_dv(const unsigned int order,
                                  const Real vel_U,
                                  const Real gamma,
                                  const Real mach);

        
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
                                            MAST::BoundaryConditionBase& bc) = 0;
        
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
                                            MAST::BoundaryConditionBase& bc) = 0;

        
        
        /*!
         *    Calculates the force vector and Jacobian due to piston-theory
         *    based surface pressure on the entire element domain.
         *    This is applicable for only 1D and 2D elements. The order
         *    of the boundary condition and direction of fluid flow are
         *    obtained from the BoundaryConditionBase object.
         */
        virtual bool piston_theory_residual_sensitivity(const MAST::FunctionBase& p,
                                                        bool request_jacobian,
                                                        RealVectorX &f,
                                                        RealMatrixX& jac_xdot,
                                                        RealMatrixX& jac,
                                                        MAST::BoundaryConditionBase& bc) = 0;
        
        /*!
         *    Calculates the force vector and Jacobian due to piston-theory
         *    based surface pressure on the side identified by \p side.
         *    The order of the boundary condition and direction of fluid
         *     flow are obtained from the BoundaryConditionBase object.
         */
        virtual bool piston_theory_residual_sensitivity(const MAST::FunctionBase& p,
                                                        bool request_jacobian,
                                                        RealVectorX &f,
                                                        RealMatrixX& jac_xdot,
                                                        RealMatrixX& jac,
                                                        const unsigned int side,
                                                        MAST::BoundaryConditionBase& bc) = 0;

        
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
         MAST::BoundaryConditionBase& bc) = 0;
        
        
        /*!
         *    Calculates the force vector and Jacobian due to surface pressure
         *    applied on the entire element domain. This is applicable for
         *    only 1D and 2D elements. The implementation can be used as 
         *    the evaluation of \f$ df(x_s, x_f)/dx_f dx_f \f$, or the 
         *    contribution of the off-diagonal Jacobian times fluid solution 
         *    perturbation.
         */
        virtual bool
        linearized_frequency_domain_surface_pressure_residual
        (bool request_jacobian,
         ComplexVectorX& f,
         ComplexMatrixX& jac,
         MAST::BoundaryConditionBase& bc);
        
        
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
         MAST::BoundaryConditionBase& bc) = 0;
        
        
        /*!
         *    Calculates the sensitivity of force vector and Jacobian due
         *    to surface pressure applied on the entire element domain. This
         *    is applicable for only 1D and 2D elements.
         */
        virtual bool
        linearized_frequency_domain_surface_pressure_residual_sensitivity
        (const MAST::FunctionBase& p,
         bool request_jacobian,
         ComplexVectorX& f,
         ComplexMatrixX& jac,
         MAST::BoundaryConditionBase& bc) {
            
            libmesh_error(); // to be implemented
        }
        
        
        /*!
         *    Calculates the force vector and Jacobian due to surface pressure.
         */
        virtual bool
        surface_pressure_residual_sensitivity(const MAST::FunctionBase& p,
                                              bool request_jacobian,
                                              RealVectorX& f,
                                              RealMatrixX& jac,
                                              const unsigned int side,
                                              MAST::BoundaryConditionBase& bc) = 0;
        
        /*!
         *    Calculates the force vector and Jacobian due to surface pressure.
         *    this should be implemented for each element type
         */
        virtual bool
        surface_pressure_residual_sensitivity(const MAST::FunctionBase& p,
                                              bool request_jacobian,
                                              RealVectorX& f,
                                              RealMatrixX& jac,
                                              MAST::BoundaryConditionBase& bc);

        
        
        /*!
         *    Calculates the force vector and Jacobian due to thermal stresses.
         *    this should be implemented for each element type
         */
        virtual bool thermal_residual(bool request_jacobian,
                                      RealVectorX& f,
                                      RealMatrixX& jac,
                                      MAST::BoundaryConditionBase& bc) = 0;
        
        /*!
         *    Calculates the sensitivity of force vector and Jacobian due to
         *    thermal stresses. this should be implemented for each element type
         */
        virtual bool thermal_residual_sensitivity(const MAST::FunctionBase& p,
                                                  bool request_jacobian,
                                                  RealVectorX& f,
                                                  RealMatrixX& jac,
                                                  MAST::BoundaryConditionBase& bc) = 0;
        
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
                                                        RealMatrixX& jac) = 0;
        
        /*!
         *   element property
         */
        const MAST::ElementPropertyCardBase& _property;
        
        
        /*!
         *   local solution
         */
        RealVectorX _local_sol;


        /*!
         *   local perturbed solution
         */
        RealVectorX _local_delta_sol;

        
        /*!
         *   local solution sensitivity
         */
        RealVectorX _local_sol_sens;


        /*!
         *   local perturbed solution sensitivity
         */
        RealVectorX _local_delta_sol_sens;

        
        /*!
         *   local velocity
         */
        RealVectorX _local_vel;

        
        /*!
         *   local perturbed velocity
         */
        RealVectorX _local_delta_vel;

        
        /*!
         *   local velocity sensitivity
         */
        RealVectorX _local_vel_sens;

        
        /*!
         *   local perturbed velocity  sensitivity
         */
        RealVectorX _local_delta_vel_sens;

        
        /*!
         *   local acceleration
         */
        RealVectorX _local_accel;

        
        /*!
         *   local perturbed acceleration
         */
        RealVectorX _local_delta_accel;

        
        /*!
         *   local acceleration sensitivity
         */
        RealVectorX _local_accel_sens;

        
        /*!
         *   local perturbed acceleration sensitivity
         */
        RealVectorX _local_delta_accel_sens;

        
        /*!
         *   incompatible mode solution vector
         */
        RealVectorX* _incompatible_sol;
        
    };
    
    
    /*!
     *    builds the structural element for the specified element type
     */
    std::unique_ptr<MAST::StructuralElementBase>
    build_structural_element(MAST::SystemInitialization& sys,
                             MAST::AssemblyBase& assembly,
                             const MAST::GeomElem& elem,
                             const MAST::ElementPropertyCardBase& p);
    
}


#endif // __mast__structural_element_base__
