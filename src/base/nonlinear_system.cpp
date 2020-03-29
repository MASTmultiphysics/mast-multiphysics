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

// C++ includes
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

// MAST includes
#include "base/nonlinear_system.h"
#include "base/physics_discipline_base.h"
#include "base/nonlinear_implicit_assembly.h"
#include "base/eigenproblem_assembly.h"
#include "base/parameter.h"
#include "base/output_assembly_elem_operations.h"
#include "solver/slepc_eigen_solver.h"

// libMesh includes
#include "libmesh/numeric_vector.h"
#include "libmesh/equation_systems.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dof_map.h"
#include "libmesh/nonlinear_solver.h"
#include "libmesh/petsc_linear_solver.h"
#include "libmesh/xdr_cxx.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/utility.h"
#include "libmesh/libmesh_version.h"
#include "libmesh/generic_projector.h"
#include "libmesh/wrapped_functor.h"
#include "libmesh/fem_context.h"
#include "libmesh/parallel.h"


MAST::NonlinearSystem::NonlinearSystem(libMesh::EquationSystems& es,
                                       const std::string& name,
                                       const unsigned int number):
libMesh::NonlinearImplicitSystem(es, name, number),
_initialize_B_matrix                  (false),
matrix_A                              (nullptr),
matrix_B                              (nullptr),
eigen_solver                          (nullptr),
_condensed_dofs_initialized           (false),
_exchange_A_and_B                     (false),
_n_requested_eigenpairs               (0),
_n_converged_eigenpairs               (0),
_n_iterations                         (0),
_is_generalized_eigenproblem          (false),
_eigen_problem_type                   (libMesh::NHEP),
_operation                            (MAST::NonlinearSystem::NONE) {
    
}



MAST::NonlinearSystem::~NonlinearSystem() {
    
    this->clear();
}





void
MAST::NonlinearSystem::clear() {
    
    
    // delete the matricies
    if (matrix_A) delete matrix_A;
    if (matrix_B) delete matrix_B;
    
    // nullptr-out the matricies.
    matrix_A = nullptr;
    matrix_B = nullptr;
    
    // clear the solver
    if (eigen_solver.get()) {
      eigen_solver->clear();
    }
    
    libMesh::NonlinearImplicitSystem::clear();
}


void
MAST::NonlinearSystem::set_eigenproblem_type (libMesh::EigenProblemType ept) {
    
    _eigen_problem_type = ept;
    
}




void
MAST::NonlinearSystem::init_data () {
    
    // initialize parent data
    libMesh::NonlinearImplicitSystem::init_data();
    
    // define the type of eigenproblem
    if (_eigen_problem_type == libMesh::GNHEP ||
        _eigen_problem_type == libMesh::GHEP  ||
        _eigen_problem_type == libMesh::GHIEP)
        _is_generalized_eigenproblem = true;
    
    
    libMesh::DofMap& dof_map = this->get_dof_map();
    
    // build the system matrix
    matrix_A = libMesh::SparseMatrix<Real>::build(this->comm()).release();
    dof_map.attach_matrix(*matrix_A);
    matrix_A->init();
    matrix_A->zero();
    
    if (_is_generalized_eigenproblem || _initialize_B_matrix) {
        
        matrix_B = libMesh::SparseMatrix<Real>::build(this->comm()).release();
        dof_map.attach_matrix(*matrix_B);
        matrix_B->init();
        matrix_B->zero();
    }
    
    eigen_solver.reset(new MAST::SlepcEigenSolver(this->comm()));
    if (libMesh::on_command_line("--solver_system_names")) {
        
        EPS eps =  eigen_solver.get()->eps();
        std::string nm = this->name() + "_";
        EPSSetOptionsPrefix(eps, nm.c_str());
    }
    eigen_solver->set_eigenproblem_type(_eigen_problem_type);
    
    
    linear_solver.reset(new libMesh::PetscLinearSolver<Real>(this->comm()));
    if (libMesh::on_command_line("--solver_system_names")) {
        
        std::string nm = this->name() + "_";
        linear_solver->init(nm.c_str());
    }
}




void MAST::NonlinearSystem::reinit () {
    // initialize parent data
    libMesh::NonlinearImplicitSystem::reinit();
    
    // Clear the matrices
    matrix_A->clear();
    
    if (_is_generalized_eigenproblem || _initialize_B_matrix)
        matrix_B->clear();
    
    eigen_solver.reset(new MAST::SlepcEigenSolver(this->comm()));
    if (libMesh::on_command_line("--solver_system_names")) {
        
        EPS eps =  eigen_solver.get()->eps();
        std::string nm = this->name() + "_";
        EPSSetOptionsPrefix(eps, nm.c_str());
    }
    eigen_solver->set_eigenproblem_type(_eigen_problem_type);
    
    
    linear_solver.reset(new libMesh::PetscLinearSolver<Real>(this->comm()));
    if (libMesh::on_command_line("--solver_system_names")) {
        
        std::string nm = this->name() + "_";
        linear_solver->init(nm.c_str());
    }

    matrix_A->init();
    matrix_A->zero();
    
    if (_is_generalized_eigenproblem || _initialize_B_matrix)
    {
        matrix_B->init();
        matrix_B->zero();
    }
}

std::pair<unsigned int, Real>
MAST::NonlinearSystem::get_linear_solve_parameters() {
    
    this->set_solver_parameters();
    return libMesh::NonlinearImplicitSystem::get_linear_solve_parameters();
}


void
MAST::NonlinearSystem::solve(MAST::AssemblyElemOperations& elem_ops,
                             MAST::AssemblyBase&  assembly) {
    
    libmesh_assert(_operation == MAST::NonlinearSystem::NONE);
    
    _operation = MAST::NonlinearSystem::NONLINEAR_SOLVE;
    assembly.set_elem_operation_object(elem_ops);
    
    libMesh::NonlinearImplicitSystem::ComputeResidualandJacobian
    *old_ptr = this->nonlinear_solver->residual_and_jacobian_object;
    
    this->nonlinear_solver->residual_and_jacobian_object = &assembly;
    //if (assembly.get_solver_monitor())
    //    assembly.get_solver_monitor()->init(assembly);
    
    libMesh::NonlinearImplicitSystem::solve();
    
    // enforce constraints on the solution since NonlinearImplicitSystem only
    // enforces the constraints on current_local_solution.
    this->get_dof_map().enforce_constraints_exactly(*this, this->solution.get());
    
    this->nonlinear_solver->residual_and_jacobian_object = old_ptr;
    
    assembly.clear_elem_operation_object();
    _operation = MAST::NonlinearSystem::NONE;
}



void
MAST::NonlinearSystem::eigenproblem_solve(MAST::AssemblyElemOperations& elem_ops,
                                          MAST::EigenproblemAssembly& assembly) {
    
    libmesh_assert(_operation == MAST::NonlinearSystem::NONE);
    
    _operation = MAST::NonlinearSystem::EIGENPROBLEM_SOLVE;

    START_LOG("eigensolve()", "NonlinearSystem");
    
    assembly.set_elem_operation_object(elem_ops);
    
    // A reference to the EquationSystems
    libMesh::EquationSystems& es = this->get_equation_systems();
    
    // check that necessary parameters have been set
    libmesh_assert(_n_requested_eigenpairs);
    
    es.parameters.set<unsigned int>("eigenpairs")    =   _n_requested_eigenpairs;
    es.parameters.set<unsigned int>("basis vectors") = 5*_n_requested_eigenpairs;
    
    // Get the tolerance for the solver and the maximum
    // number of iterations. Here, we simply adopt the linear solver
    // specific parameters.
    const Real
    tol    = es.parameters.get<Real>("linear solver tolerance");
    
    const unsigned int
    maxits = es.parameters.get<unsigned int>("linear solver maximum iterations"),
    nev    = es.parameters.get<unsigned int>("eigenpairs"),
    ncv    = es.parameters.get<unsigned int>("basis vectors");
    
    std::pair<unsigned int, unsigned int> solve_data;
    
    libMesh::SparseMatrix<Real>
    *eig_A  = nullptr,
    *eig_B  = nullptr;
    
    // assemble the matrices
    assembly.eigenproblem_assemble(matrix_A, matrix_B);

    // If we haven't initialized any condensed dofs,
    // just use the default eigen_system
    if (!_condensed_dofs_initialized) {
        
        // call the solver depending on the type of eigenproblem
        if (generalized())  {
            
            //in case of a generalized eigenproblem
            
            // exchange the matrices if requested by the user
            if (!_exchange_A_and_B) {
                eig_A  =  matrix_A;
                eig_B  =  matrix_B;
            }
            else {
                eig_B  =  matrix_A;
                eig_A  =  matrix_B;
            }
            
            solve_data = eigen_solver->solve_generalized (*eig_A,
                                                          *eig_B,
                                                          nev,
                                                          ncv,
                                                          tol,
                                                          maxits);
        }
        else {
            
            libmesh_assert (!matrix_B);
            
            //in case of a standard eigenproblem
            solve_data = eigen_solver->solve_standard (*matrix_A,
                                                       nev,
                                                       ncv,
                                                       tol,
                                                       maxits);
        }
    }
    else {
        
        // If we reach here, then there should be some non-condensed dofs
        libmesh_assert(!_local_non_condensed_dofs_vector.empty());

        std::unique_ptr<libMesh::SparseMatrix<Real> >
        condensed_matrix_A(libMesh::SparseMatrix<Real>::build(this->comm()).release()),
        condensed_matrix_B(libMesh::SparseMatrix<Real>::build(this->comm()).release());

        // Now condense the matrices
        matrix_A->create_submatrix(*condensed_matrix_A,
                                   _local_non_condensed_dofs_vector,
                                   _local_non_condensed_dofs_vector);
        
        
        if (generalized()) {
            
            matrix_B->create_submatrix(*condensed_matrix_B,
                                       _local_non_condensed_dofs_vector,
                                       _local_non_condensed_dofs_vector);
        }
        
        // call the solver depending on the type of eigenproblem
        if ( generalized() ) {
            
            //in case of a generalized eigenproblem
            
            // exchange the matrices if requested by the user
            if (!_exchange_A_and_B) {
                eig_A  =  condensed_matrix_A.get();
                eig_B  =  condensed_matrix_B.get();
            }
            else {
                eig_B  =  condensed_matrix_A.get();
                eig_A  =  condensed_matrix_B.get();
            }
            
            solve_data = eigen_solver->solve_generalized(*eig_A,
                                                         *eig_B,
                                                         nev,
                                                         ncv,
                                                         tol,
                                                         maxits);
        }
        else {
            
            libmesh_assert (!matrix_B);
            
            //in case of a standard eigenproblem
            solve_data = eigen_solver->solve_standard (*condensed_matrix_A,
                                                       nev,
                                                       ncv,
                                                       tol,
                                                       maxits);
        }
    }
    
    _n_converged_eigenpairs = solve_data.first;
    _n_iterations           = solve_data.second;
    
    assembly.clear_elem_operation_object();
    
    STOP_LOG("eigensolve()", "NonlinearSystem");
    _operation = MAST::NonlinearSystem::NONE;
}




void
MAST::NonlinearSystem::get_eigenvalue(unsigned int i, Real&  re, Real&  im) {
    
    std::pair<Real, Real>
    val   = this->eigen_solver->get_eigenvalue (i);
    
    if (!_exchange_A_and_B) {
        re   = val.first;
        im   = val.second;
    }
    else {
        Complex complex_val (val.first, val.second);
        complex_val = 1./complex_val;
        re = complex_val.real();
        im = complex_val.imag();
    }
}



void
MAST::NonlinearSystem::get_eigenpair(unsigned int i,
                                     Real&  re,
                                     Real&  im,
                                     libMesh::NumericVector<Real>& vec_re,
                                     libMesh::NumericVector<Real>* vec_im) {
    
    START_LOG("get_eigenpair()", "NonlinearSystem");
    std::pair<Real, Real>
    val;
    
    // imaginary component is needed only when the problem is non-Hermitian
    if (_eigen_problem_type == libMesh::HEP ||
        _eigen_problem_type == libMesh::GHEP)
        libmesh_assert (vec_im == nullptr);
    
    
    // If we haven't initialized any condensed dofs,
    // just use the default eigen_system
    if (!_condensed_dofs_initialized) {
        
        // call the eigen_solver get_eigenpair method
        val   = this->eigen_solver->get_eigenpair (i, vec_re, vec_im);
        
        if (!_exchange_A_and_B) {
            re   = val.first;
            im   = val.second;
        }
        else {
            Complex complex_val (val.first, val.second);
            complex_val = 1./complex_val;
            re = complex_val.real();
            im = complex_val.imag();
        }
    }
    else {
        
        // If we reach here, then there should be some non-condensed dofs
        libmesh_assert(!_local_non_condensed_dofs_vector.empty());
        
        std::unique_ptr< libMesh::NumericVector<Real> >
        temp_re(libMesh::NumericVector<Real>::build(this->comm()).release()),
        temp_im;
        
        unsigned int
        n_local   = (unsigned int)_local_non_condensed_dofs_vector.size(),
        n         = n_local;
        this->comm().sum(n);
        
        // initialize the vectors
        temp_re->init (n, n_local, false, libMesh::PARALLEL);
        
        // imaginary only if the problem is non-Hermitian
        if (vec_im) {
            
            temp_im.reset(libMesh::NumericVector<Real>::build(this->comm()).release());
            temp_im->init (n, n_local, false, libMesh::PARALLEL);
        }
        

        // call the eigen_solver get_eigenpair method
        val   = this->eigen_solver->get_eigenpair (i, *temp_re, temp_im.get());
        
        if (!_exchange_A_and_B) {
            re   = val.first;
            im   = val.second;
        }
        else {
            Complex complex_val (val.first, val.second);
            complex_val = 1./complex_val;
            re = complex_val.real();
            im = complex_val.imag();
        }

        
        // Now map temp to solution. Loop over local entries of local_non_condensed_dofs_vector
        // the real part
        vec_re.zero();
        for (unsigned int j=0; j<_local_non_condensed_dofs_vector.size(); j++) {
            
            unsigned int index = _local_non_condensed_dofs_vector[j];
            vec_re.set(index,(*temp_re)(temp_re->first_local_index()+j));
        }
        vec_re.close();
        
        // now the imaginary part if it was provided
        if (vec_im) {
            
            vec_im->zero();
            
            for (unsigned int j=0; j<_local_non_condensed_dofs_vector.size(); j++) {
                
                unsigned int index = _local_non_condensed_dofs_vector[j];
                vec_im->set(index,(*temp_im)(temp_im->first_local_index()+j));
            }
            
            vec_im->close();
        }
    }
    
    // scale the eigenvector so that it has a unit inner product
    // with the B matrix.
    switch (_eigen_problem_type) {
        case libMesh::HEP: {
            Real
            v = vec_re.dot(vec_re);
            
            // make sure that v is a nonzero value
            libmesh_assert_greater(v, 0.);
            vec_re.scale(1./std::sqrt(v));
        }
            break;
            
        case libMesh::GHEP: {
            
            std::unique_ptr<libMesh::NumericVector<Real> >
            tmp(vec_re.zero_clone().release());
            
            // inner product with respect to B matrix
            matrix_B->vector_mult(*tmp, vec_re);

            Real
            v = tmp->dot(vec_re);
            
            // make sure that v is a nonzero value
            libmesh_assert_greater(v, 0.);
            vec_re.scale(1./std::sqrt(v));
        }
            break;
            
            
        case libMesh::NHEP:   // to be implemented
        case libMesh::GNHEP:  // to be implemented
        case libMesh::GHIEP:  // to be implemented
        default:
            libmesh_error(); // should not get here
    }

    this->update();
    
    STOP_LOG("get_eigenpair()", "NonlinearSystem");
}





void
MAST::NonlinearSystem::
eigenproblem_sensitivity_solve (MAST::AssemblyElemOperations&    elem_ops,
                                MAST::EigenproblemAssembly&      assembly,
                                const MAST::FunctionBase&        f,
                                std::vector<Real>&               sens,
                                const std::vector<unsigned int>* indices) {
    
    // make sure that eigensolution is already available
    libmesh_assert(_n_converged_eigenpairs);

    assembly.set_elem_operation_object(elem_ops);
    
    // the sensitivity is calculated based on the inner product of the left and
    // right eigen vectors.
    //
    //        y^T [A] x - lambda y^T [B] x = 0
    //    where y and x are the left and right eigenvectors, respectively.
    //    Therefore,
    //        d lambda/dp = (y^T (d[A]/dp - lambda d[B]/dp) x) / (y^T [B] x)
    //
    //    the denominator remain constant for all sensitivity calculations.
    //
    const unsigned int
    nconv  = std::min(_n_requested_eigenpairs, _n_converged_eigenpairs),
    n_calc = indices?(unsigned int)indices->size():nconv;

    libmesh_assert_equal_to(n_calc, nconv);
    
    std::vector<unsigned int> indices_to_calculate;
    if (indices) {
        indices_to_calculate = *indices;
        for (unsigned int i=0; i<n_calc; i++) libmesh_assert_less(indices_to_calculate[i], nconv);
    }
    else {
        // calculate all
        indices_to_calculate.resize(n_calc);
        for (unsigned int i=0; i<n_calc; i++) indices_to_calculate[i] = i;
    }
    
    std::vector<Real>
    denom(n_calc, 0.);
    sens.resize(n_calc, 0.);
    
    std::vector<libMesh::NumericVector<Real>*>
    x_right (n_calc);
    //x_left  (_n_converged_eigenpairs);
    
    std::unique_ptr<libMesh::NumericVector<Real> >
    tmp     (this->solution->zero_clone().release());

    std::vector<Real>
    eig (n_calc);
    
    Real
    re  = 0.,
    im  = 0.;
    
    
    for (unsigned int i=0; i<n_calc; i++) {

        x_right[i] = (this->solution->zero_clone().release());
        
        
        switch (_eigen_problem_type) {
                
            case libMesh::HEP: {
                // right and left eigenvectors are same
                // imaginary part of eigenvector for real matrices is zero
                this->get_eigenpair(indices_to_calculate[i], re, im, *x_right[i], nullptr);
                denom[i] = x_right[i]->dot(*x_right[i]);               // x^H x
                eig[i]   = re;
            }
                break;
                
            case libMesh::GHEP: {
                // imaginary part of eigenvector for real matrices is zero
                this->get_eigenpair(indices_to_calculate[i], re, im, *x_right[i], nullptr);
                matrix_B->vector_mult(*tmp, *x_right[i]);
                denom[i] = x_right[i]->dot(*tmp);                  // x^H B x
                eig[i]   = re;
            }
                break;
                
            default:
                // to be implemented for the non-Hermitian problems
                libmesh_error();
                break;
        }
    }
    
    // calculate sensitivity of matrix quantities
    assembly.eigenproblem_sensitivity_assemble(f, matrix_A, matrix_B);
    
    
    // now calculate sensitivity of each eigenvalue for the parameter
    for (unsigned int i=0; i<n_calc; i++) {
        
        switch (_eigen_problem_type) {
                
            case libMesh::HEP: {
                
                matrix_A->vector_mult(*tmp, *x_right[i]);
                sens[i] = x_right[i]->dot(*tmp);                  // x^H A' x
                sens[i]-= eig[i] * x_right[i]->dot(*x_right[i]);  // - lambda x^H x
                sens[i] /= denom[i];                              // x^H x
            }
                break;
                
            case libMesh::GHEP: {
                
                matrix_A->vector_mult(*tmp, *x_right[i]);
                sens[i] = x_right[i]->dot(*tmp);              // x^H A' x
                matrix_B->vector_mult(*tmp, *x_right[i]);
                sens[i]-= eig[i] * x_right[i]->dot(*tmp);     // - lambda x^H B' x
                sens[i] /= denom[i];                          // x^H B x
            }
                break;
                
            default:
                // to be implemented for the non-Hermitian problems
                libmesh_error();
                break;
        }
    }
    
    // now delete the x_right vectors
    for (unsigned int i=0; i<x_right.size(); i++)
        delete x_right[i];
    
    assembly.clear_elem_operation_object();
}



void
MAST::NonlinearSystem::
initialize_condensed_dofs(MAST::PhysicsDisciplineBase& physics) {
    
    // This has been adapted from
    // libMesh::CondensedEigenSystem::initialize_condensed_dofs()
    const libMesh::DofMap & dof_map = this->get_dof_map();

    std::set<unsigned int> global_dirichlet_dofs_set;
    
    physics.get_system_dirichlet_bc_dofs(*this, global_dirichlet_dofs_set);
    
    
    // First, put all local dofs into non_dirichlet_dofs_set and
    std::set<unsigned int> local_non_condensed_dofs_set;
    for(unsigned int i=this->get_dof_map().first_dof();
        i<this->get_dof_map().end_dof();
        i++)
        if (!dof_map.is_constrained_dof(i))
            local_non_condensed_dofs_set.insert(i);
    
    // Now erase the condensed dofs
    std::set<unsigned int>::iterator
    iter     = global_dirichlet_dofs_set.begin(),
    iter_end = global_dirichlet_dofs_set.end();
    
    for ( ; iter != iter_end ; ++iter) {
        
        unsigned int condensed_dof_index = *iter;
        
        if ( (this->get_dof_map().first_dof() <= condensed_dof_index) &&
            (condensed_dof_index < this->get_dof_map().end_dof()) )
            local_non_condensed_dofs_set.erase(condensed_dof_index);
    }
    
    // Finally, move local_non_condensed_dofs_set over to a vector for
    // convenience in solve()
    iter     = local_non_condensed_dofs_set.begin();
    iter_end = local_non_condensed_dofs_set.end();
    
    _local_non_condensed_dofs_vector.clear();
    
    for ( ; iter != iter_end; ++iter)
        _local_non_condensed_dofs_vector.push_back(*iter);
    
    _condensed_dofs_initialized = true;
}





void
MAST::NonlinearSystem::sensitivity_solve(MAST::AssemblyElemOperations& elem_ops,
                                         MAST::AssemblyBase&           assembly,
                                         const MAST::FunctionBase&     p,
                                         bool                          if_assemble_jacobian) {

    libmesh_assert(_operation == MAST::NonlinearSystem::NONE);
    
    _operation = MAST::NonlinearSystem::FORWARD_SENSITIVITY_SOLVE;
    
    // Log how long the linear solve takes.
    LOG_SCOPE("sensitivity_solve()", "NonlinearSystem");

    assembly.set_elem_operation_object(elem_ops);
    
    libMesh::NumericVector<Real>
    &dsol  = this->add_sensitivity_solution(),
    &rhs   = this->add_sensitivity_rhs();

    if (if_assemble_jacobian)
        assembly.residual_and_jacobian(*solution, nullptr, matrix, *this);
    assembly.sensitivity_assemble(p, rhs);
    
    rhs.scale(-1.);
    
    // The sensitivity problem is linear
    // Our iteration counts and residuals will be sums of the individual
    // results
    std::pair<unsigned int, Real>
    solver_params = this->get_linear_solve_parameters();
    
    // Solve the linear system.
    libMesh::SparseMatrix<Real> * pc = this->request_matrix("Preconditioner");
    
    std::pair<unsigned int, Real> rval =
    this->linear_solver->solve (*matrix, pc,
                                dsol,
                                rhs,
                                solver_params.second,
                                solver_params.first);
    
    // The linear solver may not have fit our constraints exactly
#ifdef LIBMESH_ENABLE_CONSTRAINTS
    this->get_dof_map().enforce_constraints_exactly (*this, &dsol, /* homogeneous = */ true);
#endif
    
    assembly.clear_elem_operation_object();
    
    _operation = MAST::NonlinearSystem::NONE;
    

}



void
MAST::NonlinearSystem::adjoint_solve(MAST::AssemblyElemOperations&       elem_ops,
                                     MAST::OutputAssemblyElemOperations& output,
                                     MAST::AssemblyBase&                 assembly,
                                     bool if_assemble_jacobian) {
    

    libmesh_assert(_operation == MAST::NonlinearSystem::NONE);
    
    _operation = MAST::NonlinearSystem::ADJOINT_SOLVE;
    
    // Log how long the linear solve takes.
    LOG_SCOPE("adjoint_solve()", "NonlinearSystem");
    
    libMesh::NumericVector<Real>
    &dsol  = this->add_adjoint_solution(),
    &rhs   = this->add_adjoint_rhs();

    assembly.set_elem_operation_object(elem_ops);

    if (if_assemble_jacobian)
        assembly.residual_and_jacobian(*solution, nullptr, matrix, *this);
    
    assembly.calculate_output_derivative(*solution, output, rhs);

    assembly.clear_elem_operation_object();

    rhs.scale(-1.);
    
    // Our iteration counts and residuals will be sums of the individual
    // results
    std::pair<unsigned int, Real>
    solver_params = this->get_linear_solve_parameters();
    
    const std::pair<unsigned int, Real> rval =
    linear_solver->adjoint_solve (*matrix,
                                  dsol,
                                  rhs,
                                  solver_params.second,
                                  solver_params.first);
    
    // The linear solver may not have fit our constraints exactly
#ifdef LIBMESH_ENABLE_CONSTRAINTS
    this->get_dof_map().enforce_adjoint_constraints_exactly(dsol, 0);
#endif
    
    _operation = MAST::NonlinearSystem::NONE;
}





void
MAST::NonlinearSystem::write_out_vector(libMesh::NumericVector<Real>& vec,
                                        const std::string & directory_name,
                                        const std::string & data_name,
                                        const bool write_binary_vectors)
{
    LOG_SCOPE("write_out_vector()", "NonlinearSystem");
    
    if (this->processor_id() == 0)
    {
        // Make a directory to store all the data files
        libMesh::Utility::mkdir(directory_name.c_str());
    }
    
    // Make sure processors are synced up before we begin
    this->comm().barrier();
    
    std::ostringstream file_name;
    const std::string suffix = (write_binary_vectors ? ".xdr" : ".dat");
    
    file_name << directory_name << "/" << data_name << "_data" << suffix;
    libMesh::Xdr bf_data(file_name.str(),
                         write_binary_vectors ? libMesh::ENCODE : libMesh::WRITE);
    
    std::string version("libMesh-" + libMesh::get_io_compatibility_version());
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
    version += " with infinite elements";
#endif
    bf_data.data(version ,"# File Format Identifier");
    
    this->write_header(bf_data, /*(unused arg)*/ version, /*write_additional_data=*/false);
    
    // Following EquationSystemsIO::write, we use a temporary numbering (node major)
    // before writing out the data
    libMesh::MeshTools::Private::globally_renumber_nodes_and_elements(this->get_mesh());
    
    // Write all vectors at once.
    {
        // Note the API wants pointers to constant vectors, hence this...
        std::vector<const libMesh::NumericVector<Real> *> bf_out(1);
        bf_out[0] = &vec;
        this->write_serialized_vectors (bf_data, bf_out);
    }
    
    
    // set the current version
    bf_data.set_version(LIBMESH_VERSION_ID(LIBMESH_MAJOR_VERSION,
                                           LIBMESH_MINOR_VERSION,
                                           LIBMESH_MICRO_VERSION));
    
    
    // Undo the temporary renumbering
    this->get_mesh().fix_broken_node_and_element_numbering();
}



void
MAST::NonlinearSystem::read_in_vector(libMesh::NumericVector<Real>& vec,
                                      const std::string & directory_name,
                                      const std::string & data_name,
                                      const bool read_binary_vector) {
    
    LOG_SCOPE("read_in_vector()", "NonlinearSystem");
    
    // Make sure processors are synced up before we begin
    this->comm().barrier();
    
    
    // Following EquationSystemsIO::read, we use a temporary numbering (node major)
    // before writing out the data. For the sake of efficiency, we do this once for
    // all the vectors that we read in.
    libMesh::MeshTools::Private::globally_renumber_nodes_and_elements(this->get_mesh());
    
    std::ostringstream file_name;
    const std::string suffix = (read_binary_vector ? ".xdr" : ".dat");
    file_name.str("");
    file_name << directory_name
    << "/" << data_name
    << "_data" << suffix;
    
    // On processor zero check to be sure the file exists
    if (this->processor_id() == 0)
    {
        struct stat stat_info;
        int stat_result = stat(file_name.str().c_str(), &stat_info);
        
        if (stat_result != 0)
            libmesh_error_msg("File does not exist: " + file_name.str());
    }
    
    if (!std::ifstream(file_name.str()))
        libmesh_error_msg("File missing: " + file_name.str());

    libMesh::Xdr vector_data(file_name.str(),
                             read_binary_vector ? libMesh::DECODE : libMesh::READ);
    
    // Read the header data. This block of code is based on EquationSystems::_read_impl.
    {
        std::string version;
        vector_data.data(version);
        
        const std::string libMesh_label = "libMesh-";
        std::string::size_type lm_pos = version.find(libMesh_label);
        if (lm_pos==std::string::npos)
        {
            libmesh_error_msg("version info missing in Xdr header");
        }
        
        std::istringstream iss(version.substr(lm_pos + libMesh_label.size()));
        int ver_major = 0, ver_minor = 0, ver_patch = 0;
        char dot;
        iss >> ver_major >> dot >> ver_minor >> dot >> ver_patch;
        vector_data.set_version(LIBMESH_VERSION_ID(ver_major, ver_minor, ver_patch));
        
        // Actually read the header data. When we do this, set read_header=false
        // so taht we do not reinit sys, since we assume that it has already been
        // set up properly (e.g. the appropriate variables have already been added).
        this->read_header(vector_data, version, /*read_header=*/false, /*read_additional_data=*/false);
    }
    
    std::vector<libMesh::NumericVector<Real> *> bf_in(1);
    bf_in[0] = &vec;
    this->read_serialized_vectors (vector_data, bf_in);
    
    // Undo the temporary renumbering
    this->get_mesh().fix_broken_node_and_element_numbering();
}



void
MAST::NonlinearSystem::
project_vector_without_dirichlet (libMesh::NumericVector<Real> & new_vector,
                                  libMesh::FunctionBase<Real>& f) const {
    
    LOG_SCOPE ("project_vector_without_dirichlet()", "NonlinearSystem");
    
    libMesh::ConstElemRange active_local_range
    (this->get_mesh().active_local_elements_begin(),
     this->get_mesh().active_local_elements_end() );
    
    libMesh::VectorSetAction<Real> setter(new_vector);
    
    const unsigned int n_variables = this->n_vars();
    
    std::vector<unsigned int> vars(n_variables);
    for (unsigned int i=0; i != n_variables; ++i)
        vars[i] = i;
    
    // Use a typedef to make the calling sequence for parallel_for() a bit more readable
    typedef
    libMesh::GenericProjector<libMesh::FEMFunctionWrapper<Real>, libMesh::FEMFunctionWrapper<libMesh::Gradient>,
    Real, libMesh::VectorSetAction<Real>> FEMProjector;
    
    libMesh::WrappedFunctor<Real>     f_fem(f);
    libMesh::FEMFunctionWrapper<Real> fw(f_fem);
    
#if (LIBMESH_MAJOR_VERSION == 1 && LIBMESH_MINOR_VERSION < 5)
    libMesh::Threads::parallel_for
    (active_local_range,
     FEMProjector(*this, fw, nullptr, setter, vars));
#else
    FEMProjector projector(*this, fw, nullptr, setter, vars);
    projector.project(active_local_range);
#endif

    // Also, load values into the SCALAR dofs
    // Note: We assume that all SCALAR dofs are on the
    // processor with highest ID
    if (this->processor_id() == (this->n_processors()-1))
    {
        // FIXME: Do we want to first check for SCALAR vars before building this? [PB]
        libMesh::FEMContext context( *this );
        
        const libMesh::DofMap & dof_map = this->get_dof_map();
        for (unsigned int var=0; var<this->n_vars(); var++)
            if (this->variable(var).type().family == libMesh::SCALAR)
            {
                // FIXME: We reinit with an arbitrary element in case the user
                //        doesn't override FEMFunctionBase::component. Is there
                //        any use case we're missing? [PB]
                libMesh::Elem * el = const_cast<libMesh::Elem *>(*(this->get_mesh().active_local_elements_begin()));
                context.pre_fe_reinit(*this, el);
                
                std::vector<libMesh::dof_id_type> SCALAR_indices;
                dof_map.SCALAR_dof_indices (SCALAR_indices, var);
                const unsigned int n_SCALAR_dofs =
                libMesh::cast_int<unsigned int>(SCALAR_indices.size());
                
                for (unsigned int i=0; i<n_SCALAR_dofs; i++)
                {
                    const libMesh::dof_id_type global_index = SCALAR_indices[i];
                    const unsigned int component_index =
                    this->variable_scalar_number(var,i);
                    
                    new_vector.set(global_index, f_fem.component(context, component_index, libMesh::Point(), this->time));
                }
            }
    }
    
    new_vector.close();
}
