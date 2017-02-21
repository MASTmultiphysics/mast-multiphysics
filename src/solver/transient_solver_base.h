/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2017  Manav Bhatia
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


#ifndef __mast__transient_solver_base__
#define __mast__transient_solver_base__

// MAST includes
#include "base/mast_data_types.h"


// libMesh includes
#include "libmesh/numeric_vector.h"


namespace MAST {
    
    // Forward declerations
    class TransientAssembly;
    class ElementBase;
    class NonlinearSystem;
    
    
    class TransientSolverBase {
    public:
        TransientSolverBase();
        
        virtual ~TransientSolverBase();

        /*!
         *   Attaches the assembly object that provides the x_dot, M and J
         *   quantities for the element
         */
        void set_assembly(MAST::TransientAssembly& assembly);

        /*!
         *   Clears the assembly object
         */
        void clear_assembly();
        
        /*!
         *   time step
         */
        Real dt;

        /*!
         *    @returns the highest order time derivative that the solver 
         *    will handle
         */
        virtual int ode_order() const = 0;
        
        /*!
         *    @returns a reference to the localized solution from
         *    iteration number = current - prev_iter. So, \par prev_iter = 0
         *    gives the current solution estimate. Note that \par prev_iter
         *    cannot be greater than the total number of iterations that this
         *    solver stores solutions for.
         */
        libMesh::NumericVector<Real>&
        solution(unsigned int prev_iter = 0) const;
        
        /*!
         *    @returns a reference to the localized velocity from
         *    iteration number = current - prev_iter. So, \par prev_iter = 0
         *    gives the current velocity estimate. Note that \par prev_iter
         *    cannot be greater than the total number of iterations that this
         *    solver stores solutions for.
         */
        libMesh::NumericVector<Real>&
        velocity(unsigned int prev_iter = 0) const;

        
        /*!
         *    @returns a reference to the localized acceleration from
         *    iteration number = current - prev_iter. So, \par prev_iter = 0
         *    gives the current acceleration estimate. Note that \par prev_iter
         *    cannot be greater than the total number of iterations that this
         *    solver stores solutions for.
         */
        libMesh::NumericVector<Real>&
        acceleration(unsigned int prev_iter = 0) const;

        
        /*!
         *   solves the current time step for solution and velocity
         */
        virtual void solve() = 0;
        
        
        /*!
         *    To be used only for initial conditions.
         *    Initializes the highest derivative solution using the solution 
         *    and its low order time derivatives at the specified time step.
         *    Then advances the time step so that the solver is ready for 
         *    time integration.
         */
        void solve_highest_derivative_and_advance_time_step();


        
        /*!
         *   advances the time step and copies the current solution to old
         *   solution, and so on.
         */
        virtual void advance_time_step();

        
        /*!
         *    localizes the relevant solutions for system assembly. The
         *    calling function has to delete the pointers to these vectors
         */
        virtual void
        build_local_quantities(const libMesh::NumericVector<Real>& current_sol,
                               std::vector<libMesh::NumericVector<Real>*>& qtys);

        /*!
         *    localizes the relevant perturbations in solutions for system
         *    assembly. 
         *    \param current_sol  the perturbation in current displacement 
         *    \f$ \Delta X \f$
         *    \param qtys upon returning, this vector is populated such that the
         *    ith element of the vector is \f$ (d (d^iX/dt^i)/ dX) dX \Delta X$
         */
        virtual void
        build_perturbed_local_quantities
        (const libMesh::NumericVector<Real>& current_sol,
         std::vector<libMesh::NumericVector<Real>*>& qtys);

        /*!
         *    TransientAssembly needs to be able to call the assembly routines
         *    of this class.
         */
        friend class MAST::TransientAssembly;

        
    protected:
        
        /*!
         *    flag to check if this is the first time step.
         */
        bool  _first_step;
        
        /*!
         *    @returns the number of iterations for which solution and velocity
         *    are to be stored.
         */
        virtual unsigned int _n_iters_to_store() const = 0;
        
        /*!
         *    provides the element with the transient data for calculations
         */
        virtual void
        _set_element_data(const std::vector<libMesh::dof_id_type>& dof_indices,
                          const std::vector<libMesh::NumericVector<Real>*>& sols,
                          MAST::ElementBase& elem) = 0;

        /*!
         *    provides the element with the transient data for calculations
         */
        virtual void
        _set_element_perturbed_data
        (const std::vector<libMesh::dof_id_type>& dof_indices,
         const std::vector<libMesh::NumericVector<Real>*>& sols,
         MAST::ElementBase& elem) = 0;

        
        /*!
         *    update the transient velocity based on the current solution
         */
        virtual void _update_velocity(libMesh::NumericVector<Real>& vel,
                                      const libMesh::NumericVector<Real>& sol) = 0;
        
        /*!
         *    update the transient acceleration based on the current solution
         */
        virtual void _update_acceleration(libMesh::NumericVector<Real>& acc,
                                          const libMesh::NumericVector<Real>& sol) = 0;

        
        /*!
         *    update the perturbation in transient velocity based on the 
         *    current perturbed solution
         */
        virtual void
        _update_delta_velocity(libMesh::NumericVector<Real>& vel,
                               const libMesh::NumericVector<Real>& sol) = 0;
        
        /*!
         *    update the perturbation in transient acceleration based on the
         *    current perturbed solution
         */
        virtual void
        _update_delta_acceleration(libMesh::NumericVector<Real>& acc,
                                   const libMesh::NumericVector<Real>& sol) = 0;

        
        /*!
         *   performs the element calculations over \par elem, and returns
         *   the element vector and matrix quantities in \par mat and
         *   \par vec, respectively. \par if_jac tells the method to also
         *   assemble the Jacobian, in addition to the residual vector.
         */
        virtual void
        _elem_calculations(MAST::ElementBase& elem,
                           const std::vector<libMesh::dof_id_type>& dof_indices,
                           bool if_jac,
                           RealVectorX& vec,
                           RealMatrixX& mat) = 0;
        
        /*!
         *   performs the element calculations over \par elem, and returns
         *   the element vector quantity in \par vec. The vector quantity only
         *   include the \f$ [J] \{dX\} f$ components, so the inherited classes
         *   must ensure that no component of constant forces (traction/body
         *   forces/etc.) are added to this vector.
         */
        virtual void
        _elem_linearized_jacobian_solution_product(MAST::ElementBase& elem,
                                                   const std::vector<libMesh::dof_id_type>& dof_indices,
                                                   RealVectorX& vec) = 0;

        
        /*!
         *   performs the element sensitivity calculations over \par elem,
         *   and returns the element residual sensitivity in \par vec .
         */
        virtual void
        _elem_sensitivity_calculations(MAST::ElementBase& elem,
                                       const std::vector<libMesh::dof_id_type>& dof_indices,
                                       RealVectorX& vec) = 0;
        
        /*!
         *   Associated TransientAssembly object that provides the 
         *   element level quantities
         */
        MAST::TransientAssembly* _assembly;
        
        /*!
         *   NonlinearImplicitSystem for which this object is
         *   calculating the solution
         */
        MAST::NonlinearSystem* _system;
        
        
        /*!
         *    flag if the current procedure is to evaluate the highest ode
         *    derivative solution, or to evaluate solution at current time step.
         */
        bool   _if_highest_derivative_solution;

    };

}


#endif // __mast__transient_solver_base__
