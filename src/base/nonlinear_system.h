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

#ifndef __mast__nonlinear_system_h__
#define __mast__nonlinear_system_h__

// C++ includes
#include <memory>

// MAST includes
#include "base/mast_data_types.h"

// libMesh includes
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/enum_eigen_solver_type.h"
#include "libmesh/eigen_system.h"


namespace MAST {
    
    
    // Forward declerations
    class Parameter;
    class SlepcEigenSolver;
    class EigenproblemAssembly;
    class PhysicsDisciplineBase;
    class OutputAssemblyBase;
    
    
    /*!
     *   This class implements a system for solution of nonlinear systems. 
     *   The class also provides a mechanism to solve associated eigenproblems
     *   \f$ {\bf A} {\bf x} = \lambda {\bf B} {\bf x} \f$. The user needs to 
     *   provide an eigenproblem assembly class to assemble the A and B
     *   matrices.
     */
    class NonlinearSystem:
    public libMesh::NonlinearImplicitSystem {
        
    public:
        
        /*!
         *    Default constructor
         */
        NonlinearSystem(libMesh::EquationSystems& es,
                         const std::string& name,
                         const unsigned int number);
        
        
        virtual ~NonlinearSystem();
        
        
        /*!
         *    flag to also initialize the B matrix. Must be called before
         *    EquationsSystems::init(). This is false by default.
         */
        void set_init_B_matrix() {
            _initialize_B_matrix = true;
        }
        
        
        /*!
         * Clear all the data structures associated with
         * the system.
         */
        virtual void clear () libmesh_override;
        
        
        /**
         * Reinitializes the member data fields associated with
         * the system, so that, e.g., \p assemble() may be used.
         */
        virtual void reinit () libmesh_override;
        
        
        /**
         *   calculates and stores the sensitivity RHS in System::add_sensitivity_rhs(0).
         *   This assumes that \p parameters is of size 1. For solving sensitivity
         *   wrt multiple parameters, the user should call this method multiple 
         *   times with the parameter vector modified.
         *   The RHS is \f$ - \frac{\partial R(U,p)}{\partial p}\f$
         */
        virtual void
        assemble_residual_derivatives
        (const libMesh::ParameterVector & parameters) override;
        
        
        /**
         * Calls user qoi derivative function.
         */
        virtual void
        assemble_qoi_derivative (const libMesh::QoISet & qoi_indices = libMesh::QoISet(),
                                 bool include_liftfunc = true,
                                 bool apply_constraints = true) libmesh_override;

        
        /*!
         *   solves the adjoint problem for the provided output function
         */
        void adjoint_solve(MAST::OutputAssemblyBase& output);
        
        
        /**
         * Assembles & solves the eigen system.
         */
        virtual void eigenproblem_solve () ;
        
        /**
         * Solves the sensitivity system, for the provided parameters. The return
         * parameters are irrelevant for EigenSystem. Sensitivity of eigenvalues
         * are returned in \p sens.
         *
         * This method is only implemented in some derived classes.
         */
        virtual void
        eigenproblem_sensitivity_solve (const libMesh::ParameterVector& parameters,
                                        std::vector<Real>& sens) ;
        
        
        /*!
         *  Assembles the matrix_A and matrix_B for eigensolution
         */
        virtual void assemble_eigensystem();

        
        /*!
         *  Assembles the sensitivity of matrix_A and matrix_B with respect to the
         *  specified parameter
         */
        virtual void
        assemble_eigensystem_sensitivity(const libMesh::ParameterVector& parameters,
                                         const unsigned int p);
        
        
        /*!
         * @returns real and imaginary parts of the ith eigenvalue for the
         * eigenproblem \f$ {\bf A} {\bf x} = \lambda {\bf B} {\bf x} \f$, and
         * the associated eigenvector.
         *
         * The returned eigenvector will be scaled such that it has a unit 
         * inner product with respect to the B matrix.
         *
         * Note that eigen problem type HEP or GHEP, \p vec_im must be nullptr,
         * and for eigenproblem type NHEP or GNHEP, the real and imag
         * parts of the eigenvector are copied to \p vec_re and \p vec_im,
         * respectively. If \p vec_im is not provided, then only the real part
         * will be copied to \p vec_re.
         */
        virtual void
        get_eigenpair (unsigned int i,
                       Real&  re,
                       Real&  im,
                       libMesh::NumericVector<Real>& vec_re,
                       libMesh::NumericVector<Real>* vec_im = nullptr);
        
        /*!
         * sets the flag to exchange the A and B matrices for a generalized 
         * eigenvalue problem. This is needed typically when the B matrix is
         * not positive semi-definite.
         */
        void set_exchange_A_and_B (bool flag) {_exchange_A_and_B = flag;}

        /**
         * sets the number of eigenvalues requested
         */
        void set_n_requested_eigenvalues (unsigned int n)
        { _n_requested_eigenpairs = n; };

        
        /**
         * @returns the number of converged eigenpairs.
         */
        unsigned int
        get_n_converged_eigenvalues () const {return _n_converged_eigenpairs;}

        /**
         * @returns the number of requested eigenpairs.
         */
        unsigned int
        get_n_requested_eigenvalues () const {return _n_requested_eigenpairs;}

        /**
         * @returns the number of eigen solver iterations.
         */
        unsigned int get_n_iterations () const {return _n_iterations;}
        
        /**
         * Sets the type of the current eigen problem.
         */
        void set_eigenproblem_type (libMesh::EigenProblemType ept);
        
        /**
         * @returns the eigen problem type.
         */
        libMesh::EigenProblemType
        get_eigenproblem_type () const {return _eigen_problem_type;}
        
        /**
         * @returns true if the underlying problem is generalized
         * , false otherwise.
         */
        bool generalized () const { return _is_generalized_eigenproblem; }
        
        
        /**
         * Register a user object to use in assembling the system matrix sensitivities
         */
        void
        attach_eigenproblem_assemble_object(MAST::EigenproblemAssembly& assemble);
        
        
        /**
         *    clears the user specified function/object for use in assembling
         *    the system matrix sensitivities
         */
        void reset_eigenproblem_assemble_object();
        
        /**
         * The system matrix for standard eigenvalue problems.
         */
        libMesh::SparseMatrix<Real> *matrix_A;
        
        /**
         * A second system matrix for generalized eigenvalue problems.
         */
        libMesh::SparseMatrix<Real> *matrix_B;
        
        
        /**
         * The EigenSolver, definig which interface, i.e solver
         * package to use.
         */
        std::unique_ptr<MAST::SlepcEigenSolver> eigen_solver;
        
        /**
         * Loop over the dofs on each processor to initialize the list
         * of non-condensed dofs. These are the dofs in the system that
         * are not contained in \p global_dirichlet_dofs_set.
         */
        void initialize_condensed_dofs(MAST::PhysicsDisciplineBase& physics);
        
        /**
         * @return the global number of non-condensed dofs in the system.
         */
        unsigned int n_global_non_condensed_dofs() const;
        
    protected:
        
        
        /**
         * Initializes the member data fields associated with
         * the system, so that, e.g., \p assemble() may be used.
         */
        virtual void init_data () libmesh_override;
        
        
        /**
         * Set the _n_converged_eigenpairs member, useful for
         * subclasses of EigenSystem.
         */
        void set_n_converged (unsigned int nconv)
        { _n_converged_eigenpairs = nconv; }
        
        /**
         * Set the _n_iterations member, useful for subclasses of
         * EigenSystem.
         */
        void set_n_iterations (unsigned int its)
        { _n_iterations = its;}
        
        
        /*!
         *   initialize the B matrix in addition to A, which might be needed
         *   for solution of complex system of equations using PC field split
         */
        bool _initialize_B_matrix;
        
        
        /**
         * A private flag to indicate whether the condensed dofs
         * have been initialized.
         */
        bool                               _condensed_dofs_initialized;
        
        /**
         * The number of requested eigenpairs.
         */
        unsigned int                       _n_requested_eigenpairs;
        
        /*!
         *  flag to exchange the A and B matrices in the eigenproblem solution
         */
        bool                               _exchange_A_and_B;
        
        /**
         * The number of converged eigenpairs.
         */
        unsigned int                       _n_converged_eigenpairs;
        
        /**
         * The number of iterations of the eigen solver algorithm.
         */
        unsigned int                       _n_iterations;
        
        /**
         * A boolean flag to indicate whether we are dealing with
         * a generalized eigenvalue problem.
         */
        bool                               _is_generalized_eigenproblem;
        
        /**
         * The type of the eigenvalue problem.
         */
        libMesh::EigenProblemType          _eigen_problem_type;
        
        /**
         * Object that assembles the sensitivity of eigen_system.
         */
        MAST::EigenproblemAssembly*        _eigenproblem_assemble_system_object;

        /*!
         *    OutputAssemblyBase object for which the adjoint calculation
         *    is being solved
         */
        MAST::OutputAssemblyBase*       _output;
        
        /**
         * Vector storing the local dof indices that will not be condensed.
         * All dofs that are not in this vector will be eliminated from
         * the system when we perform a solve.
         */
        std::vector<libMesh::dof_id_type>  _local_non_condensed_dofs_vector;
        
    };
}


#endif // __mast__nonlinear_system_h__
