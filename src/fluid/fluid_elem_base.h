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

#ifndef __mast__fluid_elem_base__
#define __mast__fluid_elem_base__

// C++ includes
#include <ostream>
#include <map>

// MAST includes
#include "base/mast_data_types.h"
#include "numerics/fem_operator_matrix.h"
//#include "BoundaryConditions/boundary_surface_motion.h"

// libMesh includes
#include "libmesh/fe_base.h"



namespace MAST {

    // Forward declerations
    class FlightCondition;
    class PrimitiveSolution;
    template <typename ValType> class SmallPerturbationPrimitiveSolution;
    
    /*!
     *   enumeration of the primitive fluid variables
     */
    enum FluidPrimitiveVars {
        RHO_PRIM,
        VEL1,
        VEL2,
        VEL3,
        TEMP
    };
    

    /*!
     *   enumeration of the conservative fluid variables
     */
    enum FluidConservativeVars {
        RHO_CONS,
        RHOVEL1,
        RHOVEL2,
        RHOVEL3,
        ETOT
    };
    
    
    /*!
     *   This class provides the necessary functions to evaluate the flux 
     *   vectors and their Jacobians for both inviscid and viscous flows.
     */
    class FluidElemBase {
        
    public:
        
        /*!
         *   Constructor
         */
        FluidElemBase():
        _if_viscous(false),
        _include_pressure_switch(false),
        //surface_motion(NULL),
        flight_condition(NULL),
        dim(0),
        _dissipation_scaling(1.)
        { }
        

        virtual ~FluidElemBase();
        
        
        void init_data();
        
        
        void get_infinity_vars( RealVectorX& vars_inf ) const;
        
        /*!
         *    This defines the surface motion for use with the nonlinear
         *    fluid solver. This can be used to define either a time-dependent
         *    motion, or a steady-state motion.
         */
        //MAST::SurfaceMotionBase* surface_motion;
        
        MAST::FlightCondition* flight_condition;
        
        unsigned int dim;
        
        bool if_viscous() const {
            return _if_viscous;
        }
        
        void calculate_dxidX (const unsigned int qp,
                              const libMesh::FEBase& fe,
                              RealMatrixX& dxi_dX,
                              RealMatrixX& dX_dxi);
        
        
        void
        update_solution_at_quadrature_point(const unsigned int qp,
                                            const libMesh::FEBase& fe,
                                            const RealVectorX& elem_solution,
                                            RealVectorX& conservative_sol,
                                            MAST::PrimitiveSolution& primitive_sol,
                                            MAST::FEMOperatorMatrix& B_mat,
                                            std::vector<MAST::FEMOperatorMatrix>& dB_mat);
        
        
        void
        calculate_advection_flux(const unsigned int calculate_dim,
                                 const MAST::PrimitiveSolution& sol,
                                 RealVectorX& flux);
        
        void
        calculate_diffusion_flux(const unsigned int calculate_dim,
                                 const MAST::PrimitiveSolution& sol,
                                 const RealMatrixX& stress_tensor,
                                 const RealVectorX& temp_gradient,
                                 RealVectorX& flux);
        
        /*!
         *    calculates and returns the stress tensor in \p stress_tensor.
         */
        void
        calculate_diffusion_tensors(const RealVectorX& elem_sol,
                                    const std::vector<MAST::FEMOperatorMatrix>& dB_mat,
                                    const RealMatrixX& dprim_dcons,
                                    const MAST::PrimitiveSolution& psol,
                                    RealMatrixX& stress_tensor,
                                    RealVectorX& temp_gradient);
        
        void
        calculate_conservative_variable_jacobian(const MAST::PrimitiveSolution& sol,
                                                 RealMatrixX& dcons_dprim,
                                                 RealMatrixX& dprim_dcons);
        
        void
        calculate_advection_flux_jacobian(const unsigned int calculate_dim,
                                          const MAST::PrimitiveSolution& sol,
                                          RealMatrixX& mat);
        
        void
        calculate_advection_flux_jacobian_rho_derivative(const unsigned int calculate_dim,
                                                         const MAST::PrimitiveSolution& sol,
                                                         RealMatrixX& mat);
        

        void
        calculate_advection_flux_jacobian_u1_derivative(const unsigned int calculate_dim,
                                                        const MAST::PrimitiveSolution& sol,
                                                        RealMatrixX& mat);
        
        void
        calculate_advection_flux_jacobian_u2_derivative(const unsigned int calculate_dim,
                                                        const MAST::PrimitiveSolution& sol,
                                                        RealMatrixX& mat);

        
        void
        calculate_advection_flux_jacobian_u3_derivative(const unsigned int calculate_dim,
                                                        const MAST::PrimitiveSolution& sol,
                                                        RealMatrixX& mat);
        

        void
        calculate_advection_flux_jacobian_T_derivative(const unsigned int calculate_dim,
                                                       const MAST::PrimitiveSolution& sol,
                                                       RealMatrixX& mat);
        
        
        template <typename ValType>
        void
        calculate_advection_flux_jacobian_vec_mult_second_derivative(const unsigned int calculate_dim,
                                                                     const MAST::PrimitiveSolution &sol,
                                                                     const ValType& lin_sol,
                                                                     RealMatrixX &mat);
        
        template <typename ValType>
        void
        calculate_advection_flux_jacobian_vec_adjoint_mult_second_derivative(const unsigned int calculate_dim,
                                                                             const MAST::PrimitiveSolution &sol,
                                                                             const ValType& lin_sol,
                                                                             RealMatrixX &mat);
        
        
        void calculate_diffusion_flux_jacobian(const unsigned int flux_dim,
                                               const unsigned int deriv_dim,
                                               const MAST::PrimitiveSolution& sol,
                                               RealMatrixX& mat);
        
        void calculate_advection_flux_jacobian_sensitivity_for_conservative_variable
        (const unsigned int calculate_dim,
         const MAST::PrimitiveSolution& sol,
         std::vector<RealMatrixX >& mat);
        
        
        void calculate_advection_flux_jacobian_sensitivity_for_primitive_variable
        (const unsigned int calculate_dim,
         const unsigned int primitive_var,
         const MAST::PrimitiveSolution& sol,
         RealMatrixX& mat);
        
        
        void calculate_advection_left_eigenvector_and_inverse_for_normal
        (const MAST::PrimitiveSolution& sol,
         const libMesh::Point& normal,
         RealVectorX& eig_vals,
         RealMatrixX& l_eig_mat,
         RealMatrixX& l_eig_mat_inv_tr);
        
        
        void calculate_advection_left_eigenvector_and_inverse_rho_derivative_for_normal
        (const MAST::PrimitiveSolution& sol,
         const libMesh::Point& normal,
         RealMatrixX& eig_vals,
         RealMatrixX& l_eig_mat,
         RealMatrixX& l_eig_mat_inv_tr);
        
        
        void calculate_advection_left_eigenvector_and_inverse_u1_derivative_for_normal
        (const MAST::PrimitiveSolution& sol,
         const libMesh::Point& normal,
         RealMatrixX& eig_vals,
         RealMatrixX& l_eig_mat,
         RealMatrixX& l_eig_mat_inv_tr);
        
        
        void calculate_advection_left_eigenvector_and_inverse_u2_derivative_for_normal
        (const MAST::PrimitiveSolution& sol,
         const libMesh::Point& normal,
         RealMatrixX& eig_vals,
         RealMatrixX& l_eig_mat,
         RealMatrixX& l_eig_mat_inv_tr);
        
        
        void calculate_advection_left_eigenvector_and_inverse_u3_derivative_for_normal
        (const MAST::PrimitiveSolution& sol,
         const libMesh::Point& normal,
         RealMatrixX& eig_vals,
         RealMatrixX& l_eig_mat,
         RealMatrixX& l_eig_mat_inv_tr);
        
        
        void calculate_advection_left_eigenvector_and_inverse_T_derivative_for_normal
        (const MAST::PrimitiveSolution& sol,
         const libMesh::Point& normal,
         RealMatrixX& eig_vals,
         RealMatrixX& l_eig_mat,
         RealMatrixX& l_eig_mat_inv_tr);
        
        
        template <typename ValType>
        void
        calculate_advection_outflow_flux_jacobian_vec_mult_second_derivative(const MAST::PrimitiveSolution &sol,
                                                                             const libMesh::Point &normal,
                                                                             const ValType& lin_sol,
                                                                             RealMatrixX &mat);
        
        
        template <typename ValType>
        void
        calculate_advection_outflow_flux_jacobian_vec_adjoint_mult_second_derivative(const MAST::PrimitiveSolution &sol,
                                                                                     const libMesh::Point &normal,
                                                                                     const ValType& lin_sol,
                                                                                     RealMatrixX &mat);
        
        
        void calculate_entropy_variable_jacobian(const MAST::PrimitiveSolution& sol,
                                                 RealMatrixX& dUdV,
                                                 RealMatrixX& dVdU);
        
        
        void
        calculate_advection_flux_jacobian_for_moving_solid_wall_boundary(const MAST::PrimitiveSolution& sol,
                                                                         const Real ui_ni,
                                                                         const libMesh::Point& nvec,
                                                                         const RealVectorX& dnvec,
                                                                         RealMatrixX& mat);
        
        
        void
        calculate_advection_flux_jacobian_rho_derivative_for_moving_solid_wall_boundary
        (const MAST::PrimitiveSolution& sol,
         const Real ui_ni,
         const libMesh::Point& nvec,
         const RealVectorX& dnvec,
         RealMatrixX& mat);
        
        void
        calculate_advection_flux_jacobian_u1_derivative_for_moving_solid_wall_boundary
        (const MAST::PrimitiveSolution& sol,
         const Real ui_ni,
         const libMesh::Point& nvec,
         const RealVectorX& dnvec,
         RealMatrixX& mat);
        
        void
        calculate_advection_flux_jacobian_u2_derivative_for_moving_solid_wall_boundary
        (const MAST::PrimitiveSolution& sol,
         const Real ui_ni,
         const libMesh::Point& nvec,
         const RealVectorX& dnvec,
         RealMatrixX& mat);
        
        void
        calculate_advection_flux_jacobian_u3_derivative_for_moving_solid_wall_boundary
        (const MAST::PrimitiveSolution& sol,
         const Real ui_ni,
         const libMesh::Point& nvec,
         const RealVectorX& dnvec,
         RealMatrixX& mat);
        
        void
        calculate_advection_flux_jacobian_T_derivative_for_moving_solid_wall_boundary
        (const MAST::PrimitiveSolution& sol,
         const Real ui_ni,
         const libMesh::Point& nvec,
         const RealVectorX& dnvec,
         RealMatrixX& mat);
        
        
        template <typename ValType>
        void
        calculate_advection_flux_jacobian_vec_mult_second_derivative_for_moving_solid_wall_boundary
        (const MAST::PrimitiveSolution &sol,
         const Real ui_ni,
         const libMesh::Point &nvec,
         const RealVectorX &dnvec,
         const ValType& lin_sol,
         RealMatrixX &mat);
        
        template <typename ValType>
        void
        calculate_advection_flux_jacobian_vec_adjoint_mult_second_derivative_for_moving_solid_wall_boundary
        (const MAST::PrimitiveSolution &sol,
         const Real ui_ni,
         const libMesh::Point &nvec,
         const RealVectorX &dnvec,
         const ValType& lin_sol,
         RealMatrixX &mat);
        
        
        
        void calculate_advection_flux_spatial_derivative
        (const unsigned int i, RealVectorX* flux,
         RealMatrixX* dflux_dU);
        
        
        
        void calculate_diffusive_flux_jacobian(unsigned int div_coord,
                                               unsigned int flux_coord,
                                               RealMatrixX& mat);
        
        
        bool calculate_barth_tau_matrix(const unsigned int qp,
                                        const libMesh::FEBase& fe,
                                        const MAST::PrimitiveSolution& sol,
                                        RealMatrixX& tau,
                                        std::vector<RealMatrixX >& tau_sens);
        
        bool calculate_aliabadi_tau_matrix(const unsigned int qp,
                                           const libMesh::FEBase& fe,
                                           const MAST::PrimitiveSolution& sol,
                                           RealMatrixX& tau,
                                           std::vector<RealMatrixX >& tau_sens);
        
        void calculate_hartmann_discontinuity_operator
        (const unsigned int qp,
         const libMesh::FEBase& fe,
         const MAST::PrimitiveSolution& sol,
         const RealVectorX& elem_solution,
         const std::vector<MAST::FEMOperatorMatrix>& dB_mat,
         const RealMatrixX& Ai_Bi_advection,
         RealVectorX& discontinuity_val);
        
        
        void calculate_aliabadi_discontinuity_operator
        (const unsigned int qp,
         const libMesh::FEBase& fe,
         const MAST::PrimitiveSolution& sol,
         const RealVectorX& elem_solution,
         const std::vector<MAST::FEMOperatorMatrix>& dB_mat,
         const RealMatrixX& Ai_Bi_advection,
         RealVectorX& discontinuity_val);
        
        
        template <typename ValType>
        void calculate_small_disturbance_aliabadi_discontinuity_operator
        (const unsigned int qp,
         const libMesh::FEBase& fe,
         const MAST::PrimitiveSolution& sol,
         const MAST::SmallPerturbationPrimitiveSolution<ValType>& dsol,
         const RealVectorX& elem_solution,
         const std::vector<MAST::FEMOperatorMatrix>& dB_mat,
         const RealMatrixX& Ai_Bi_advection,
         RealVectorX& discontinuity_val);
        
        
        void calculate_differential_operator_matrix
        (const unsigned int qp,
         const libMesh::FEBase& fe,
         const RealVectorX& elem_solution,
         const MAST::PrimitiveSolution& sol,
         const MAST::FEMOperatorMatrix& B_mat,
         const std::vector<MAST::FEMOperatorMatrix>& dB_mat,
         const std::vector<RealMatrixX >& Ai_advection,
         const RealMatrixX& Ai_Bi_advection,
         const std::vector<std::vector<RealMatrixX > >& Ai_sens,
         RealMatrixX& LS_operator,
         RealMatrixX& LS_sens);
        
    protected:
        
        std::vector<MAST::FluidPrimitiveVars> _active_primitive_vars;
        
        std::vector<FluidConservativeVars> _active_conservative_vars;
        
        bool _if_viscous;
        
        bool _include_pressure_switch;
        
        Real _dissipation_scaling;
    };
    
    
}

#endif /* defined(__mast__fluid_elem_base__) */
