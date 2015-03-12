//
//  structural_element_base.h
//  mast
//
//  Created by Manav Bhatia on 2/17/15.
//  Copyright (c) 2015 MAST. All rights reserved.
//

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
    
    
    
    
    class StructuralElementBase:
    public MAST::ElementBase
    {
    public:
        /*!
         *   Constructor
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
         *   side external force contribution to system residual
         */
        template <typename ValType>
        bool side_external_residual (bool request_jacobian,
                                     RealVectorX& f,
                                     RealMatrixX& jac,
                                     std::multimap<libMesh::boundary_id_type, MAST::BoundaryConditionBase*>& bc);
        
        /*!
         *   prestress force contribution to system residual
         */
        virtual bool prestress_residual (bool request_jacobian,
                                         RealVectorX& f,
                                         RealMatrixX& jac) = 0;
        
        /*!
         *   volume external force contribution to system residual
         */
        template <typename ValType>
        bool volume_external_residual (bool request_jacobian,
                                       RealVectorX& f,
                                       RealMatrixX& jac,
                                       std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*>& bc);
        
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
                                                   RealMatrixX& jac,
                                                   std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*>& bc);
        
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
         *    this should be implemented for each element type
         */
        virtual bool surface_pressure_residual(bool request_jacobian,
                                               RealVectorX& f,
                                               RealMatrixX& jac,
                                               const unsigned int side,
                                               MAST::BoundaryConditionBase& bc);
        
        
        /*!
         *    Calculates the force vector and Jacobian due to surface pressure
         *    applied on the entire element domain. This is applicable for
         *    only 1D and 2D elements.
         */
        virtual bool surface_pressure_residual(bool request_jacobian,
                                               RealVectorX& f,
                                               RealMatrixX& jac,
                                               MAST::BoundaryConditionBase& bc);
        
        
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
        virtual bool surface_pressure_residual_sensitivity(bool request_jacobian,
                                                           RealVectorX& f,
                                                           RealMatrixX& jac,
                                                           const unsigned int side,
                                                           MAST::BoundaryConditionBase& bc)
        { libmesh_error();}
        
        
        /*!
         *    Calculates the force vector and Jacobian due to surface pressure.
         *    this should be implemented for each element type
         */
        virtual bool surface_pressure_residual_sensitivity(bool request_jacobian,
                                                           RealVectorX& f,
                                                           RealMatrixX& jac,
                                                           MAST::BoundaryConditionBase& bc)
        { libmesh_error();}
        
        
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
