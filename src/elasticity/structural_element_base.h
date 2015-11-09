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

#ifndef __mast__structural_element_base__
#define __mast__structural_element_base__

// C++ includes
#include <memory>
#include <map>

// MAST includes
#include "base/elem_base.h"


namespace MAST {
    
    // Forward declerations
    class ElementPropertyCardBase;
    class LocalElemBase;
    class BoundaryConditionBase;
    class FEMOperatorMatrix;
    class OutputFunctionBase;
    
    
    class StructuralElementBase:
    public MAST::ElementBase
    {
    public:
        /*!
         *   Constructor.
         */
        StructuralElementBase(MAST::SystemInitialization& sys,
                              const libMesh::Elem& elem,
                              const MAST::ElementPropertyCardBase& p);
        
        virtual ~StructuralElementBase();
        
        
        /*!
         *   stores \p vec as solution for element level calculations,
         *   or its sensitivity if \p if_sens is true.
         */
        virtual void set_solution(const RealVectorX& vec,
                                  bool if_sens = false);
        
        
        /*!
         *    stores \p vec as velocity for element level calculations,
         *    or its sensitivity if \p if_sens is true.
         */
        virtual void set_velocity(const RealVectorX& vec,
                                  bool if_sens = false);

        
        /*!
         *    stores \p vec as acceleration for element level calculations,
         *    or its sensitivity if \p if_sens is true.
         */
        virtual void set_acceleration(const RealVectorX& vec,
                                      bool if_sens = false);

        
        /*!
         *   This is used for cases where a linearized problem is solved
         *   about a stationary base solution. This method stores
         *   \p vec as the base solution, or its sensitivity if \p
         *   if_sens is true.
         */
        virtual void set_base_solution(const RealVectorX& vec,
                                       bool if_sens = false);

        
        /*!
         *  @returns a constant reference to the element solution 
         *  (or its derivative if \par if_sens is true) in the local
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
                                        RealMatrixX& jac,
                                        bool if_ignore_ho_jac) = 0;
        
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
                                       RealMatrixX& jac)
        { libmesh_assert(false);}
        
        /*!
         *   inertial force contribution to system residual
         */
        virtual bool inertial_residual (bool request_jacobian,
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
        template <typename ValType>
        bool side_external_residual (bool request_jacobian,
                                     RealVectorX& f,
                                     RealMatrixX& jac_xdot,
                                     RealMatrixX& jac,
                                     std::multimap<libMesh::boundary_id_type, MAST::BoundaryConditionBase*>& bc);

        
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
        template <typename ValType>
        bool volume_external_residual (bool request_jacobian,
                                       RealVectorX& f,
                                       RealMatrixX& jac_xdot,
                                       RealMatrixX& jac,
                                       std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*>& bc);
        
        /*!
         *   evaluates an output quantity requested in the map over the
         *   boundary of the element that may coincide with the boundary
         *   identified in the map. The derivative with respect to the
         *   state variables is provided if \p request_derivative is true.
         */
        virtual bool volume_output_quantity (bool request_derivative,
                                             bool request_sensitivity,
                                             std::multimap<libMesh::subdomain_id_type, MAST::OutputFunctionBase*>& output);

        
        /*!
         *   sensitivity of the internal force contribution to system residual
         */
        virtual bool internal_residual_sensitivity (bool request_jacobian,
                                                    RealVectorX& f,
                                                    RealMatrixX& jac,
                                                    bool if_ignore_ho_jac) = 0;

        /*!
         *   sensitivity of the damping force contribution to system residual
         */
        virtual bool damping_residual_sensitivity (bool request_jacobian,
                                                   RealVectorX& f,
                                                   RealMatrixX& jac)
        { libmesh_assert(false);}
        
        /*!
         *   sensitivity of the inertial force contribution to system residual
         */
        virtual bool inertial_residual_sensitivity (bool request_jacobian,
                                                    RealVectorX& f,
                                                    RealMatrixX& jac)
        { libmesh_assert(false);}
        
        /*!
         *   sensitivity of the side external force contribution to system residual
         */
        template <typename ValType>
        bool side_external_residual_sensitivity (bool request_jacobian,
                                                 RealVectorX& f,
                                                 RealMatrixX& jac_xdot,
                                                 RealMatrixX& jac,
                                                 std::multimap<libMesh::boundary_id_type, MAST::BoundaryConditionBase*>& bc);
        
        /*!
         *   sensitivity of the prestress force contribution to system residual
         */
        virtual bool prestress_residual_sensitivity (bool request_jacobian,
                                                     RealVectorX& f,
                                                     RealMatrixX& jac) = 0;
        
        /*!
         *   sensitivity of the volume external force contribution to system residual
         */
        template <typename ValType>
        bool volume_external_residual_sensitivity (bool request_jacobian,
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
         *   @returns the piston theory cp value based on the 
         *   specified normal velocity and the flow parameters
         */
        Real piston_theory_cp(const unsigned int order,
                              const Real vel_normal,
                              const Real a_inf,
                              const Real gamma,
                              const Real mach);

        
        
        /*!
         *   @returns the derivative of piston theory cp value with respect 
         *   to the normal velocity. 
         */
        Real piston_theory_dcp_dvn(const unsigned int order,
                                   const Real vel_normal,
                                   const Real a_inf,
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
        virtual bool piston_theory_residual_sensitivity(bool request_jacobian,
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
        virtual bool piston_theory_residual_sensitivity(bool request_jacobian,
                                                        RealVectorX &f,
                                                        RealMatrixX& jac_xdot,
                                                        RealMatrixX& jac,
                                                        const unsigned int side,
                                                        MAST::BoundaryConditionBase& bc) = 0;

        
        /*!
         *    Calculates the force vector and Jacobian due to small
         *    perturbation surface pressure.
         */
        template <typename ValType>
        bool
        small_disturbance_surface_pressure_residual(bool request_jacobian,
                                                    RealVectorX& f,
                                                    RealMatrixX& jac,
                                                    const unsigned int side,
                                                    MAST::BoundaryConditionBase& bc);
        
        
        /*!
         *    Calculates the force vector and Jacobian due to surface pressure
         *    applied on the entire element domain. This is applicable for
         *    only 1D and 2D elements.
         */
        template <typename ValType>
        bool
        small_disturbance_surface_pressure_residual(bool request_jacobian,
                                                    RealVectorX& f,
                                                    RealMatrixX& jac,
                                                    MAST::BoundaryConditionBase& bc);
        
        
        /*!
         *    Calculates the sensitivity of force vector and Jacobian due to small
         *     is applicable for perturbation surface pressure.
         */
        template <typename ValType>
        bool
        small_disturbance_surface_pressure_residual_sensitivity(bool request_jacobian,
                                                                RealVectorX& f,
                                                                RealMatrixX& jac,
                                                                const unsigned int side,
                                                                MAST::BoundaryConditionBase& bc);
        
        
        /*!
         *    Calculates the sensitivity of force vector and Jacobian due
         *    to surface pressure applied on the entire element domain. This
         *    is applicable for only 1D and 2D elements.
         */
        template <typename ValType>
        bool
        small_disturbance_surface_pressure_residual_sensitivity(bool request_jacobian,
                                                                RealVectorX& f,
                                                                RealMatrixX& jac,
                                                                MAST::BoundaryConditionBase& bc);
        
        
        /*!
         *    Calculates the force vector and Jacobian due to surface pressure.
         */
        virtual bool
        surface_pressure_residual_sensitivity(bool request_jacobian,
                                              RealVectorX& f,
                                              RealMatrixX& jac,
                                              const unsigned int side,
                                              MAST::BoundaryConditionBase& bc) = 0;
        
        /*!
         *    Calculates the force vector and Jacobian due to surface pressure.
         *    this should be implemented for each element type
         */
        virtual bool
        surface_pressure_residual_sensitivity(bool request_jacobian,
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
        virtual bool thermal_residual_sensitivity(bool request_jacobian,
                                                  RealVectorX& f,
                                                  RealMatrixX& jac,
                                                  MAST::BoundaryConditionBase& bc) = 0;
        
        
        /*!
         *    Calculates the stress tensor. If derivative and sensitivity
         *    with respect to the parameter \p sesitivity_param are calculated
         *    and provided if the respective flags are true.
         */
        virtual bool calculate_stress(bool request_derivative,
                                      bool request_sensitivity,
                                      MAST::OutputFunctionBase& output) = 0;

        
        
        /*!
         *   element property
         */
        const MAST::ElementPropertyCardBase& _property;
        
        
        /*!
         *   local solution
         */
        RealVectorX _local_sol;
        
        
        /*!
         *   local solution sensitivity
         */
        RealVectorX _local_sol_sens;
        
        
        /*!
         *   local velocity
         */
        RealVectorX _local_vel;
        
        
        /*!
         *   local velocity
         */
        RealVectorX _local_vel_sens;
        
        
        /*!
         *   local acceleration
         */
        RealVectorX _local_accel;
        
        
        /*!
         *   local acceleration
         */
        RealVectorX _local_accel_sens;
        
        
        /*!
         *   base solution about which a linearized solution is performed
         */
        RealVectorX _local_base_sol;
        
        
        /*!
         *   base solution sensitivity
         */
        RealVectorX _local_base_sol_sens;

        /*!
         *   incompatible mode solution vector
         */
        RealVectorX* _incompatible_sol;
        
    };
    
    
    /*!
     *    builds the structural element for the specified element type
     */
    std::auto_ptr<MAST::StructuralElementBase>
    build_structural_element(MAST::SystemInitialization& sys,
                             const libMesh::Elem& elem,
                             const MAST::ElementPropertyCardBase& p);
}


#endif // __mast__structural_element_base__
