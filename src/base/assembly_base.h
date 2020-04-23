/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2020  Manav Bhatia and MAST authors
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

#ifndef __mast__system_assembly_base_h__
#define __mast__system_assembly_base_h__

// C++ includes
#include <map>
#include <memory>


// MAST includes
#include "base/mast_data_types.h"


// libMesh includes
#include "libmesh/system.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/sparse_matrix.h"


namespace MAST {
    
    // Forward declerations
    class PhysicsDisciplineBase;
    class SystemInitialization;
    class ElementBase;
    class MeshFieldFunction;
    class NonlinearSystem;
    class FEBase;
    class AssemblyElemOperations;
    class OutputAssemblyElemOperations;
    class FunctionBase;
    
    class AssemblyBase:
    public libMesh::NonlinearImplicitSystem::ComputeResidualandJacobian {
    public:
        
        /*!
         *  constructor takes a reference to the discipline that provides
         *  the boundary conditions, volume loads, properties, etc.
         */
        AssemblyBase();
        
        /*!
         *   virtual destructor
         */
        virtual ~AssemblyBase();
        
        /*!
         *  subdomain ids for which residuakl and Jacobian contributions will not be computed. Instead,
         *  a small diagonal value will be added to the Jacobian for the dofs corresponding to this
         *  element.
         */
        std::set<unsigned int> diagonal_elem_subdomain_id;
        

        class SolverMonitor {
        public:
            SolverMonitor(){}
            virtual ~SolverMonitor(){}
            virtual void init(MAST::AssemblyBase& assembly) = 0;
            virtual void clear() = 0;
        };
        
        /*!
         *   Inherited objects from this class can be provided by the user
         *   provide assessment of whether or not an element is influenced by
         *   a give parameter. If provided, the sensitivity assembly funcitons
         *   will use this to avoid computing sensitivity on elements where the
         *   parameter has no influence. For some cases in direct sensitivity
         *   a user may want to avoid computing this data even if a solution
         *   sensitivity vector is provided when this solution is known to be
         *   zero on the element. In this case, the override flag can be
         *   set to \p true.
         */
        class ElemParameterDependence {
        public:
            ElemParameterDependence(bool o_flag): override_flag(o_flag) {}
            virtual ~ElemParameterDependence() {}
            virtual bool if_elem_depends_on_parameter(const libMesh::Elem& e,
                                                      const MAST::FunctionBase& p) const = 0;
            /*!
             *  if \p true, assume zero solution sensitivity when elem does
             *  not dependent on parameter. This can be useful for spatial
             *  association of parameter, for example, in topology optimization.
             */
            bool override_flag;
        };
        
        /*!
         *   flag to control the closing fo the Jacobian after assembly
         */
        bool close_matrix;
        
        /*!
         *   @returns a const reference to the PhysicsDisciplineBase object
         *   associated with this object
         */
        const MAST::PhysicsDisciplineBase& discipline() const;
        
        /*!
         *   @returns a non-const reference to the PhysicsDisciplineBase object
         *   associated with this object
         */
        MAST::PhysicsDisciplineBase& discipline();
        
        
        /*!
         *   @returns a reference to the element operations object
         */
        MAST::AssemblyElemOperations& get_elem_ops();
        
        
        /*!
         *   attaches a system to this discipline
         */
        virtual void
        set_discipline_and_system(MAST::PhysicsDisciplineBase& discipline,
                                  MAST::SystemInitialization& system);

        /*!
         *   This object, if provided by user, will be used to reduce
         *   unnecessary computations in sensitivity analysis assembly operations.
         *   This association is cleard when \p clear_discipline_and_system()
         *   is called.
         */
        void
        attach_elem_parameter_dependence_object(MAST::AssemblyBase::ElemParameterDependence& dep);

        
        void
        clear_elem_parameter_dependence_object();

        
        /*!
         *   clears association with a system to this discipline
         */
        virtual void
        clear_discipline_and_system();


        /*!
         *   attaches a element operation to this object, and associated
         *   \p this with the element operation object.
         */
        virtual void
        set_elem_operation_object(MAST::AssemblyElemOperations& elem_ops);

        
        /*!
         *   clears the association of \p this object with the assembly
         *   element operation object.
         */
        virtual void
        clear_elem_operation_object();

        
        /*!
         *   @returns a const reference to the MAST::NonlinearSystem object
         *   associated with this object
         */
        const MAST::NonlinearSystem& system() const;
        
        /*!
         *   @returns a non-const reference to the MAST::NonlinearSystem object
         *   associated with this object
         */
        MAST::NonlinearSystem& system();

        /*!
         *   @returns a non-const reference to the MAST::SystemInitialization
         */
        MAST::SystemInitialization& system_init();


        /*!
         *   attaches the solver monitor, which is a user provided routine
         *   that is called each time
         */
        void set_solver_monitor(MAST::AssemblyBase::SolverMonitor& monitor);

        /*!
         *   @returns the solver monitor, which is a user provided routine
         *   that is called each time
         */
        MAST::AssemblyBase::SolverMonitor* get_solver_monitor();

        /*!
         *   clears the monitor object
         */
        void clear_solver_monitor();
        
        /*!
         *   tells the assembly object that this function is will
         *   need to be initialized before each residual evaluation
         */
        void attach_solution_function(MAST::MeshFieldFunction& f);

        
        /*!
         *   removes the attachment of the solution function
         */
        void detach_solution_function();
        
        /*!
         *    function that assembles the matrices and vectors quantities for
         *    nonlinear solution
         */
        virtual void
        residual_and_jacobian (const libMesh::NumericVector<Real>& X,
                               libMesh::NumericVector<Real>* R,
                               libMesh::SparseMatrix<Real>*  J,
                               libMesh::NonlinearImplicitSystem& S) {
            libmesh_assert(false); // implement in the derived class
        }

        /**
         * Assembly function.  This function will be called
         * to assemble the RHS of the sensitivity equations (which is -1 times
         * sensitivity of system residual) prior to a solve and must
         * be provided by the user in a derived class. The method provides dR/dp
         * for \p f parameter.
         *
         * If the routine is not able to provide sensitivity for this parameter,
         * then it should return false, and the system will attempt to use
         * finite differencing.
         */
        virtual bool
        sensitivity_assemble (const libMesh::NumericVector<Real>& X,
                              bool if_localize_sol,
                              const MAST::FunctionBase& f,
                              libMesh::NumericVector<Real>& sensitivity_rhs) {
            libmesh_assert(false); // implemented in the derived class
        }

        /*!
         *   calculates the value of quantity \f$ q(X,p) \f$.
         */
        virtual void
        calculate_output(const libMesh::NumericVector<Real>& X,
                         bool if_localize_sol,
                         MAST::OutputAssemblyElemOperations& output);

        
        /*!
         *   calculates \f$ \frac{\partial q(X, p)}{\partial X} \f$
         */
        virtual void
        calculate_output_derivative(const libMesh::NumericVector<Real>& X,
                                    bool if_localize_sol,
                                    MAST::OutputAssemblyElemOperations& output,
                                    libMesh::NumericVector<Real>& dq_dX);

        
        /*!
         *   evaluates the sensitivity of the outputs in the attached
         *   discipline with respect to the parametrs in \p params.
         *   The base solution should be provided in \p X. If total sensitivity
         *   is desired, then \p dXdp should contain the sensitivity of
         *   solution wrt the parameter \p p, otherwise it can be set to
         *   nullptr. If \p dXdp is zero, the calculated sensitivity will be
         *   the partial derivarive of \p output wrt \p p.
         */
        virtual void
        calculate_output_direct_sensitivity(const libMesh::NumericVector<Real>& X,
                                            bool if_localize_sol,
                                            const libMesh::NumericVector<Real>* dXdp,
                                            bool if_localize_sol_sens,
                                            const MAST::FunctionBase& p,
                                            MAST::OutputAssemblyElemOperations& output);

        
        /*!
         *   Evaluates the total sensitivity of \p output wrt \p p using
         *   the adjoint solution provided in \p dq_dX for a linearization
         *   about solution \p X.
         */
        virtual Real
        calculate_output_adjoint_sensitivity(const libMesh::NumericVector<Real>& X,
                                             bool if_localize_sol,
                                             const libMesh::NumericVector<Real>& dq_dX,
                                             const MAST::FunctionBase& p,
                                             MAST::AssemblyElemOperations&       elem_ops,
                                             MAST::OutputAssemblyElemOperations& output,
                                             const bool include_partial_sens = true);

        
        /*!
         *   localizes the parallel vector so that the local copy
         *   stores all values necessary for calculation of the
         *   element quantities
         */
        std::unique_ptr<libMesh::NumericVector<Real> >
        build_localized_vector(const libMesh::System& sys,
                               const libMesh::NumericVector<Real>& global) const;
        
        
    protected:
        
        /*!
         *   provides assembly elem operations for use by this class
         */
        MAST::AssemblyElemOperations*  _elem_ops;

        
        /*!
         *   PhysicsDisciplineBase object for which this class is assembling
         */
        MAST::PhysicsDisciplineBase* _discipline;
        
        
        /*!
         *   System for which this assembly is performed
         */
        MAST::SystemInitialization* _system;
        
        /*!
         *   system solution that will be initialized before each solution
         */
        MAST::MeshFieldFunction* _sol_function;
        
        /*!
         *   User provided solver monitor is attached to the linear
         *   nonlinear solvers, if provided
         */
        MAST::AssemblyBase::SolverMonitor *_solver_monitor;
        
        /*!
         *   If provided by user, this object is used by sensitiivty analysis
         *   to check for whether or the current design parameter influences
         *   an element. This can be used to enhance computational efficiency.
         */
        MAST::AssemblyBase::ElemParameterDependence *_param_dependence;
    };
        
}

#endif // __mast__system_assembly_base_h__

