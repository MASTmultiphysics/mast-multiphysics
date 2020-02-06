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

#ifndef __mast__stabilized_first_order_newmark_transient_sensitivity_solver__
#define __mast__stabilized_first_order_newmark_transient_sensitivity_solver__

// MAST includes
#include "solver/transient_solver_base.h"


namespace MAST {
    
    // Forward declerations
    class OutputAssemblyElemOperations;
    
    /*!
     *   This solver implements the stabilized sensitivity analysis solver for
     *   chaotic systems where the linearized system can be unstable.
     */
    class StabilizedFirstOrderNewmarkTransientSensitivitySolver:
    public MAST::TransientSolverBase {
    public:
        StabilizedFirstOrderNewmarkTransientSensitivitySolver();
        
        virtual ~StabilizedFirstOrderNewmarkTransientSensitivitySolver();
        
        /*!
         *    \f$ \bar{a} \f$ parameter used by this solver.
         */
        Real max_amp;
        
        Real beta;
        
        /*!
         *    index of solution that is used for current linearization
         */
        unsigned int max_index;
        
        /*!
         *   sets if the eigenvalue-based stabilization will be used.
         */
        void set_eigenvalue_stabilization(bool f);
        
        /*!
         *   sets the directory where the nonlinear solutions are stored. The
         *   name of the solution is assumed to be file_root + std::string(index)
         */
        void set_nolinear_solution_location(std::string& file_root,
                                            std::string& dir);
        
        
        /*!
         *    solvers the current time step for sensitivity wrt \p f
         */
        virtual void sensitivity_solve(MAST::AssemblyBase& assembly,
                                       const MAST::FunctionBase& f);
        
        virtual Real
        evaluate_q_sens_for_previous_interval(MAST::AssemblyBase& assembly,
                                              const MAST::FunctionBase& p,
                                              MAST::OutputAssemblyElemOperations& output);
        
        virtual void
        set_element_data(const std::vector<libMesh::dof_id_type>& dof_indices,
                         const std::vector<libMesh::NumericVector<Real>*>& sols);
        
        
        virtual void
        extract_element_sensitivity_data(const std::vector<libMesh::dof_id_type>& dof_indices,
                                         const std::vector<libMesh::NumericVector<Real>*>& sols,
                                         std::vector<RealVectorX>& local_sols);
        
        virtual void
        update_sensitivity_velocity(libMesh::NumericVector<Real>&       vec,
                                    const libMesh::NumericVector<Real>& sol);
        
        virtual void
        elem_calculations(bool if_jac,
                          RealVectorX& vec,
                          RealMatrixX& mat);
        
        virtual void
        elem_sensitivity_calculations(const MAST::FunctionBase& f,
                                      RealVectorX& vec);
        
        virtual void
        elem_sensitivity_contribution_previous_timestep(const std::vector<RealVectorX>& prev_sols,
                                                        RealVectorX& vec);
        
        virtual void
        update_velocity(libMesh::NumericVector<Real>& vel,
                        const libMesh::NumericVector<Real>& sol);
        
        virtual void
        update_acceleration(libMesh::NumericVector<Real>& acc,
                                         const libMesh::NumericVector<Real>& sol){
            libmesh_error();
        }
        
        virtual void
        update_sensitivity_acceleration(libMesh::NumericVector<Real>& acc,
                                                     const libMesh::NumericVector<Real>& sol){
            libmesh_error();
        }
        
        virtual void
        update_delta_velocity(libMesh::NumericVector<Real>& vel,
                              const libMesh::NumericVector<Real>& sol){
            libmesh_error();
        }
        
        virtual void
        update_delta_acceleration(libMesh::NumericVector<Real>& acc,
                                  const libMesh::NumericVector<Real>& sol){
            libmesh_error();
        }

        virtual void
        elem_linearized_jacobian_solution_product(RealVectorX& vec) {
            libmesh_error();
        }

        virtual void
        elem_shape_sensitivity_calculations(const MAST::FunctionBase& f,
                                            RealVectorX& vec)  {
            libmesh_error();
        }

        virtual void
        elem_topology_sensitivity_calculations(const MAST::FunctionBase& f,
                                               RealVectorX& vec) {
            libmesh_error();
        }


        virtual void
        elem_topology_sensitivity_calculations(const MAST::FunctionBase& f,
                                               const MAST::FieldFunction<RealVectorX>& vel,
                                               RealVectorX& vec) {
            libmesh_error();
        }

        virtual void
        elem_second_derivative_dot_solution_assembly(RealMatrixX& mat)  {
            libmesh_error();
        }

        virtual void
        set_element_perturbed_data
        (const std::vector<libMesh::dof_id_type>& dof_indices,
         const std::vector<libMesh::NumericVector<Real>*>& sols) {
            libmesh_error();
        }
    protected:

        Real _compute_norm_amplification_factor(const libMesh::NumericVector<Real>& sol0,
                                                const libMesh::NumericVector<Real>& sol1);
        
        Real _compute_eig_amplification_factor(libMesh::SparseMatrix<Real>& A,
                                               libMesh::SparseMatrix<Real>& B);


        bool         _use_eigenvalue_stabilization;
        bool         _assemble_mass;
        Real         _t0;
        unsigned int _index0, _index1;
        
        // nonlinear solutions are stored in this directory
        std::string  _sol_name_root, _sol_dir;
    };
    
}

#endif // __mast__stabilized_first_order_newmark_transient_sensitivity_solver__

