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
        
        
        class SolverMonitor {
        public:
            SolverMonitor(){}
            virtual ~SolverMonitor(){}
            virtual void init(MAST::AssemblyBase& assembly) = 0;
            virtual void clear() = 0;
        };
        
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
         *   attaches a system to this discipline, and vice-a-versa
         */
        virtual void
        attach_discipline_and_system(MAST::AssemblyElemOperations& elem_ops,
                                     MAST::PhysicsDisciplineBase& discipline,
                                     MAST::SystemInitialization& system);

        
        /*!
         *   Reattaches to the same system that was attached earlier.
         *
         *   This cannot be called if the clear_discipline_and_system() method
         *   has been called.
         */
        virtual void
        reattach_to_system() = 0;

        
        /*!
         *   clears association with a system to this discipline, and vice-a-versa
         */
        virtual void
        clear_discipline_and_system( );

        
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
         * for \par f parameter.
         *
         * If the routine is not able to provide sensitivity for this parameter,
         * then it should return false, and the system will attempt to use
         * finite differencing.
         */
        virtual bool
        sensitivity_assemble (const MAST::FunctionBase& f,
                              libMesh::NumericVector<Real>& sensitivity_rhs) {
            libmesh_assert(false); // implemented in the derived class
        }

        /*!
         *   calculates the value of quantity \f$ q(X,p) \f$.
         */
        virtual void
        calculate_output(const libMesh::NumericVector<Real>& X,
                         MAST::OutputAssemblyElemOperations& output);

        
        /*!
         *   calculates \f$ \frac{\partial q(X, p)}{\partial X} \f$
         */
        virtual void
        calculate_output_derivative(const libMesh::NumericVector<Real>& X,
                                    MAST::OutputAssemblyElemOperations& output,
                                    libMesh::NumericVector<Real>& dq_dX);

        
        /*!
         *   evaluates the sensitivity of the outputs in the attached
         *   discipline with respect to the parametrs in \par params.
         *   The base solution should be provided in \par X. If total sensitivity
         *   is desired, then \par dXdp should contain the sensitivity of
         *   solution wrt the parameter \par p. If this \par dXdp is zero,
         *   the calculated sensitivity will be the partial derivarive of
         *   \par output wrt \par p.
         */
        virtual void
        calculate_output_direct_sensitivity(const libMesh::NumericVector<Real>& X,
                                            const libMesh::NumericVector<Real>& dXdp,
                                            const MAST::FunctionBase& p,
                                            MAST::OutputAssemblyElemOperations& output);

        
        /*!
         *   Evaluates the total sensitivity of \par output wrt \par p using
         *   the adjoint solution provided in \par dq_dX for a linearization
         *   about solution \par X.
         */
        virtual void
        calculate_output_adjoint_sensitivity(const libMesh::NumericVector<Real>& X,
                                             const libMesh::NumericVector<Real>& dq_dX,
                                             const MAST::FunctionBase& p,
                                             MAST::OutputAssemblyElemOperations& output);

        
        /*!
         *   localizes the parallel vector so that the local copy
         *   stores all values necessary for calculation of the
         *   element quantities
         */
        std::unique_ptr<libMesh::NumericVector<Real> >
        build_localized_vector(const libMesh::System& sys,
                               const libMesh::NumericVector<Real>& global) const;
        
        
        /*!
         *   @returns a MAST::FEBase object for calculation of finite element
         *   quantities. For all standard applications this is a wrapper
         *   around the libMesh::FEBase class, which is specialized for
         *   cut-cell applications where a sub-finite element is created
         *   for element integration.
         */
        virtual std::unique_ptr<MAST::FEBase>
        build_fe(const libMesh::Elem& e);
        

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
    };
        
}

#endif // __mast__system_assembly_base_h__

