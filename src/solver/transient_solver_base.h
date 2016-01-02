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


#ifndef __mast__transient_solver_base__
#define __mast__transient_solver_base__

// MAST includes
#include "base/mast_data_types.h"

// libMesh includes
#include "libmesh/nonlinear_implicit_system.h"


namespace MAST {
    
    // Forward declerations
    class TransientAssembly;
    class ElementBase;
    
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
         *    @returns a const reference to the localized solution from 
         *    iteration number = current - prev_iter. So, \par prev_iter = 0
         *    gives the current solution estimate. Note that \par prev_iter
         *    cannot be greater than the total number of iterations that this
         *    solver stores solutions for.
         */
        const libMesh::NumericVector<Real>&
        solution(unsigned int prev_iter = 0) const;
        
        /*!
         *    @returns a const reference to the localized velocity from
         *    iteration number = current - prev_iter. So, \par prev_iter = 0
         *    gives the current velocity estimate. Note that \par prev_iter
         *    cannot be greater than the total number of iterations that this
         *    solver stores solutions for.
         */
        const libMesh::NumericVector<Real>&
        velocity(unsigned int prev_iter = 0) const;

        
        /*!
         *    @returns a const reference to the localized acceleration from
         *    iteration number = current - prev_iter. So, \par prev_iter = 0
         *    gives the current acceleration estimate. Note that \par prev_iter
         *    cannot be greater than the total number of iterations that this
         *    solver stores solutions for.
         */
        const libMesh::NumericVector<Real>&
        acceleration(unsigned int prev_iter = 0) const;

        
        /*!
         *   solves the current time step for solution and velocity
         */
        virtual void solve() = 0;
        
        
        /*!
         *   This should be called before the first time step. All state vectors
         *   from the current and the previous time steps are set equal to
         *   \par val. This is useful for system with scalar unknowns.
         */
        virtual void set_initial_condition(Real val);

        
        /*!
         *   advances the time step and copies the current solution to old
         *   solution, and so on.
         */
        virtual void advance_time_step();

        
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
         *    localizes the relevant solutions for system assembly.
         */
        virtual void _localize_solutions(const libMesh::NumericVector<Real>& current_sol);
        
        /*!
         *    provides the element with the transient data for calculations
         */
        virtual void _set_element_data(std::vector<libMesh::dof_id_type>& dof_indices,
                                       MAST::ElementBase& elem) = 0;
        
        /*!
         *    update the transient velocity based on the current solution
         */
        virtual void _update_velocity(libMesh::NumericVector<Real>& vec) = 0;
        
        /*!
         *    update the transient acceleration based on the current solution
         */
        virtual void _update_acceleration(libMesh::NumericVector<Real>& vec) = 0;
        
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
        libMesh::NonlinearImplicitSystem* _system;
        
        /*!
         *   localized solution vector from previous time step needed for
         *   element calculations on local processor
         */
        std::vector<libMesh::NumericVector<Real>*> _solution;
        
        
        /*!
         *   localized velocity vector needed for element calculations
         *   on local processor
         */
        std::vector<libMesh::NumericVector<Real>*> _velocity;

        /*!
         *   localized acceleration vector needed for element calculations
         *   on local processor
         */
        std::vector<libMesh::NumericVector<Real>*> _acceleration;

    };

}


#endif // __mast__transient_solver_base__
