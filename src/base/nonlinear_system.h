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
    class AssemblyBase;
    class PhysicsDisciplineBase;
    class AssemblyElemOperations;
    class OutputAssemblyElemOperations;
    class FunctionBase;
    class EigenproblemAssembly;
    
    
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

        
        enum Operation {
            NONLINEAR_SOLVE,
            EIGENPROBLEM_SOLVE,
            FORWARD_SENSITIVITY_SOLVE,
            ADJOINT_SOLVE,
            NONE
        };
        
        
        /*!
         *   @returns the current operation of the system
         */
        MAST::NonlinearSystem::Operation operation() {
            
            return _operation;
        }

        
        /*!
         *   sets the current operation of the system
         */
        void set_operation(MAST::NonlinearSystem::Operation op) {
            
            _operation = op;
        }

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


        /*!
         *   calls NonlinearImplicitSystem::set_solver_parameters() before
         *   accessing the values.
         */
        virtual std::pair<unsigned int, Real>
        get_linear_solve_parameters();

        
        /*!
         *  solves the nonlinear problem with the specified assembly operation
         *  object
         */
        virtual void solve(MAST::AssemblyElemOperations& elem_ops,
                           MAST::AssemblyBase&           assembly);
        
        
        /*!
         *   Solves the sensitivity problem for the provided parameter.
         *   The Jacobian will be assembled before adjoint solve if
         *   \p if_assemble_jacobian is \p true.
         */
        virtual void sensitivity_solve(const libMesh::NumericVector<Real>& X,
                                       bool if_localize_sol,
                                       MAST::AssemblyElemOperations&   elem_ops,
                                       MAST::AssemblyBase&             assembly,
                                       const MAST::FunctionBase&       p,
                                       bool if_assemble_jacobian = true);

        
        /*!
         *   solves the adjoint problem for the provided output function.
         *   The Jacobian will be assembled before adjoint solve if
         *   \p if_assemble_jacobian is \p true.
         */
        virtual void adjoint_solve(const libMesh::NumericVector<Real>& X,
                                   bool if_localize_sol,
                                   MAST::AssemblyElemOperations&       elem_ops,
                                   MAST::OutputAssemblyElemOperations& output,
                                   MAST::AssemblyBase&                 assembly,
                                   bool if_assemble_jacobian           = true);
        
        
        /**
         * Assembles & solves the eigen system.
         */
        virtual void eigenproblem_solve(MAST::AssemblyElemOperations& elem_ops,
                                        MAST::EigenproblemAssembly&   assembly);
        
        /**
         * Solves the sensitivity system, for the provided parameters.
         * Sensitivity of eigenvalues are returned in \p sens. This is more
         * If only a subset of sensitivities are needed, then the indices
         * can be passed in the last argument. If the last argument is not
         * provided, then sensitivity of all eigenvalues will be computed and
         * returned in \p sens. 
         */
        virtual void
        eigenproblem_sensitivity_solve (MAST::AssemblyElemOperations& elem_ops,
                                        MAST::EigenproblemAssembly& assembly,
                                        const MAST::FunctionBase& f,
                                        std::vector<Real>& sens,
                                        const std::vector<unsigned int>* indices=nullptr);

        
        /*!
         * gets the real and imaginary parts of the ith eigenvalue for the
         * eigenproblem \f$ {\bf A} {\bf x} = \lambda {\bf B} {\bf x} \f$, and
         * the associated eigenvector.
         */
        virtual void
        get_eigenvalue (unsigned int i, Real&  re, Real&  im);

        
        /*!
         * gets the real and imaginary parts of the ith eigenvalue for the
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
        
        /*!
         *  The LinearSolver for solution of the linear equations
         */
        std::unique_ptr<libMesh::LinearSolver<Real>> linear_solver;
        
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
        
        
        /*!
         *   writes the specified vector with the specified name in a directory.
         */
        void write_out_vector(libMesh::NumericVector<Real>& vec,
                              const std::string & directory_name,
                              const std::string & data_name,
                              const bool write_binary_vectors);

        
        /*!
         *   reads the specified vector with the specified name in a directory.
         */
        void read_in_vector(libMesh::NumericVector<Real>& vec,
                            const std::string & directory_name,
                            const std::string & data_name,
                            const bool read_binary_vectors);

        void
        project_vector_without_dirichlet (libMesh::NumericVector<Real> & new_vector,
                                          libMesh::FunctionBase<Real>& f) const;
        
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
        
        /*!
         *   current operation of the system
         */
        MAST::NonlinearSystem::Operation  _operation;
        
        /**
         * Vector storing the local dof indices that will not be condensed.
         * All dofs that are not in this vector will be eliminated from
         * the system when we perform a solve.
         */
        std::vector<libMesh::dof_id_type>  _local_non_condensed_dofs_vector;
        
    };
}


#endif // __mast__nonlinear_system_h__
