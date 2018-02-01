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


#ifndef __mast__transient_solver_base__
#define __mast__transient_solver_base__

// MAST includes
#include "base/mast_data_types.h"
#include "base/nonlinear_implicit_assembly_elem_operations.h"

// libMesh includes
#include "libmesh/numeric_vector.h"


namespace MAST {
    
    // Forward declerations
    class TransientAssemblyElemOperations;
    class ElementBase;
    class NonlinearSystem;
    
    
    class TransientSolverBase:
    public MAST::NonlinearImplicitAssemblyElemOperations {
    public:
        TransientSolverBase();
        
        virtual ~TransientSolverBase();

        /*!
         *   Attaches the assembly elem operations object that provides the
         *   x_dot, M and J quantities for the element
         */
        virtual void set_assembly_ops(MAST::TransientAssemblyElemOperations& assembly_ops);

        /*!
         *   Clears the assembly elem operations object
         */
        virtual void clear_assembly_ops();
        
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
         *   @returns a reference to the localized solution sensitivity. The
         *   \p prev_iter parameter is similar to the solution() method.
         */
        libMesh::NumericVector<Real>&
        solution_sensitivity(unsigned int prev_iter = 0) const;

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
         *   @returns a reference to the localized velocity sensitivity. The
         *   \p prev_iter parameter is similar to the velocity() method.
         */
        libMesh::NumericVector<Real>&
        velocity_sensitivity(unsigned int prev_iter = 0) const;

        
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
         *   @returns a reference to the localized acceleration sensitivity. The
         *   \p prev_iter parameter is similar to the acceleration() method.
         */
        libMesh::NumericVector<Real>&
        acceleration_sensitivity(unsigned int prev_iter = 0) const;

        
        /*!
         *    update the transient velocity based on the current solution
         */
        virtual void update_velocity(libMesh::NumericVector<Real>& vel,
                                     const libMesh::NumericVector<Real>& sol) = 0;
        
        /*!
         *    update the transient acceleration based on the current solution
         */
        virtual void update_acceleration(libMesh::NumericVector<Real>& acc,
                                         const libMesh::NumericVector<Real>& sol) = 0;
        
        
        /*!
         *    update the perturbation in transient velocity based on the
         *    current perturbed solution
         */
        virtual void
        update_delta_velocity(libMesh::NumericVector<Real>& vel,
                              const libMesh::NumericVector<Real>& sol) = 0;
        
        /*!
         *    update the perturbation in transient acceleration based on the
         *    current perturbed solution
         */
        virtual void
        update_delta_acceleration(libMesh::NumericVector<Real>& acc,
                                  const libMesh::NumericVector<Real>& sol) = 0;

        /*!
         *   solves the current time step for solution and velocity
         */
        virtual void solve() = 0;
        
        /*!
         *    solvers the current time step for sensitivity wrt \p f
         */
        virtual void sensitivity_solve(const MAST::FunctionBase& f) = 0;

        
        /*!
         *    To be used only for initial conditions.
         *    Initializes the highest derivative solution using the solution 
         *    and its low order time derivatives at the specified time step.
         *    Then advances the time step so that the solver is ready for 
         *    time integration.
         */
        void solve_highest_derivative_and_advance_time_step();
        
        /*!
         *    solves for the sensitivity of highest derivative and advances
         *    the time-step.
         */
        void
        solve_highest_derivative_and_advance_time_step_with_sensitivity(const MAST::FunctionBase& f);


        
        /*!
         *   advances the time step and copies the current solution to old
         *   solution, and so on.
         */
        virtual void advance_time_step();
        
        /*!
         *   advances the time step and copies the current sensitivity solution
         *   to old sensitivity solution, and so on. Does not influence the
         *   primary solution.
         */
        virtual void advance_time_step_with_sensitivity();

        
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
         *    provides the element with the transient data for calculations
         */
        virtual void
        set_element_data(const std::vector<libMesh::dof_id_type>& dof_indices,
                         const std::vector<libMesh::NumericVector<Real>*>& sols) = 0;
        
        
        /*!
         *    provides the element with the sensitivity of transient data for
         *    calculations
         */
        virtual void
        set_element_sensitivity_data(const std::vector<libMesh::dof_id_type>& dof_indices,
                                     const std::vector<libMesh::NumericVector<Real>*>& sols) = 0;

        /*!
         *    provides the element with the transient data for calculations
         */
        virtual void
        set_element_perturbed_data
        (const std::vector<libMesh::dof_id_type>& dof_indices,
         const std::vector<libMesh::NumericVector<Real>*>& sols) = 0;

        
        /*!
         *   calls the method from TransientAssemblyElemOperations
         */
        virtual void clear_elem();
        

        virtual bool
        if_use_local_elem() const {
            libmesh_assert(false); // should not get called.
        }
        
        virtual void
        set_local_fe_data(MAST::LocalElemFE& fe,
                          const libMesh::Elem& e) const {
            libmesh_assert(false); // should not get called.
        }

    protected:
        
        /*!
         *    flag to check if this is the first time step.
         */
        bool  _first_step, _first_sensitivity_step;
        
        /*!
         *    @returns the number of iterations for which solution and velocity
         *    are to be stored.
         */
        virtual unsigned int _n_iters_to_store() const = 0;
        
        
        /*!
         *   Associated TransientAssembly object that provides the 
         *   element level quantities
         */
        MAST::TransientAssemblyElemOperations* _assembly_ops;
        
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
