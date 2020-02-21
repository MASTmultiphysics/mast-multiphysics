/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2019  Manav Bhatia
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

#ifndef __mast__warping_assembly__
#define __mast__warping_assembly__

// MAST includes
#include "base/assembly_base.h"

// libMesh includes
#include "libmesh/nonlinear_implicit_system.h"

struct geometric_properties
{
    Real A = 0.0;
    Real Qx = 0.0;
    Real Qy = 0.0;
    Real xc = 0.0;
    Real yc = 0.0;
    Real Ixxc = 0.0;
    Real Iyyc = 0.0;
    Real Ixyc = 0.0;
    Real Ixx = 0.0;
    Real Iyy = 0.0;
    Real Ixy = 0.0;
    Real Ip = 0.0;
    Real I11 = 0.0;
    Real I22 = 0.0;
    Real phi_p = 0.0;
    Real rx = 0.0;
    Real ry = 0.0;;
};


struct warping_properties
{
    Real J = 0.0;
    Real xs = 0.0;
    Real ys = 0.0;
    Real xs_t = 0.0;
    Real ys_t = 0.0;
    Real gamma = 0.0;
    Real kappa_x = 0.0;
    Real kappa_y = 0.0;
    Real kappa_xy = 0.0;
    Real Ixw = 0.0;
    Real Iyw = 0.0;
    Real Qw = 0.0;
    Real Iw = 0.0;
};

namespace MAST {
    
    // Forward declerations
    class NonlinearImplicitAssemblyElemOperations;
    
    
    class WarpingAssembly:
    public MAST::AssemblyBase {
    public:
        
        
        /*!
         *    user-provided object to perform actions
         *    after assembly and before returning to the solver. Use 
         *    \p set_post_assembly_object to provide a pointer to the object.
         */
        class PostAssemblyOperation{
            
        public:
            PostAssemblyOperation() {}
            virtual ~PostAssemblyOperation() {}
            virtual void post_assembly(const libMesh::NumericVector<Real>& X,
                                       libMesh::NumericVector<Real>* R,
                                       libMesh::SparseMatrix<Real>*  J,
                                       libMesh::NonlinearImplicitSystem& S) = 0;
        };

        /*!
         *   constructor associates this assembly object with the system
         */
        WarpingAssembly();
        
        
        /*!
         *   destructor resets the association of this assembly object with
         *   the system
         */
        virtual ~WarpingAssembly();
        
        
        
        /*!
         *   L2 norm of the last-assembled residual
         */
        Real res_l2_norm() const { return _res_l2_norm; }
        
        Real first_iter_res_l2_norm() const { return _first_iter_res_l2_norm; }

        /*!
         *   reset L2 norm of the last-assembled residual
         */
        void reset_residual_norm_history() {
            _res_l2_norm = 0.;
            _first_iter_res_l2_norm = -1.;
        }

        /*!
         *    sets the PostAssemblyOperation object for use after assembly. 
         *    Note that calling \p clear_discipline_and_system() will 
         *    clear this pointer and the user will have to call this function 
         *    again.
         */
        void
        set_post_assembly_operation(MAST::WarpingAssembly::PostAssemblyOperation& post);
        
        /*!
         *    function that assembles the matrices and vectors quantities for
         *    nonlinear solution
         */
        virtual void
        residual_and_jacobian (const libMesh::NumericVector<Real>& X,
                               libMesh::NumericVector<Real>* R,
                               libMesh::SparseMatrix<Real>*  J,
                               libMesh::NonlinearImplicitSystem& S);
        
        
        /*!
         *    calculates the product of the Jacobian and a perturbation in solution 
         *    vector \f$ [J] \{\Delta X\}  \f$. For a single discipline system the
         *    solution vector and linearized solution provided here are used. For
         *    a multiphysics system, the user must ensure that all relevant 
         *    multidisciplinary data-structures are initialized before calling 
         *    this method.
         */
        virtual void
        linearized_jacobian_solution_product(const libMesh::NumericVector<Real>& X,
                                             const libMesh::NumericVector<Real>& dX,
                                             libMesh::NumericVector<Real>& JdX,
                                             libMesh::NonlinearImplicitSystem& S);


        /*!
         *    calculates \f$ d ([J] \{\Delta X\})/ dX  \f$.
         */
        virtual void
        second_derivative_dot_solution_assembly(const libMesh::NumericVector<Real>& X,
                                                const libMesh::NumericVector<Real>& dX,
                                                libMesh::SparseMatrix<Real>& d_JdX_dX,
                                                libMesh::NonlinearImplicitSystem& S);

        
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
        sensitivity_assemble (const MAST::FunctionBase& f,
                              libMesh::NumericVector<Real>& sensitivity_rhs);
        
        /**
         * Forces the systems matrix to be symmetric using A = (A + A^T)/2
         */
        virtual void make_matrix_symmetric(libMesh::SparseMatrix<Real>* J);
        
        virtual void set_force_jacobian_symmetry(bool tf);
        
        virtual const bool get_force_jacobian_symmetry() const;
        
//         void cross_section_properties (const libMesh::NumericVector<Real>& W,
//                                        const libMesh::NumericVector<Real>& R,
//                                        const libMesh::SparseMatrix<Real>&  J);
        
        void get_loads(libMesh::NumericVector<Real>& F_warp, 
                       libMesh::NumericVector<Real>& F_shearx,
                       libMesh::NumericVector<Real>& F_sheary,
                       const Real xc, const Real yc, const Real Ixxc, 
                       const Real Iyyc, const Real Ixyc);
        
        const geometric_properties calculate_geometric_properties() const;
        
        const warping_properties calculate_warping_properties(
            const libMesh::NumericVector<Real>& F_warp, 
            const libMesh::NumericVector<Real>& Omega, 
            const libMesh::NumericVector<Real>& Psi, 
            const libMesh::NumericVector<Real>& Phi, 
            const Real A, const Real Ixxc, const Real Iyyc, const Real Ixyc,
            const Real xc, const Real yc) const;
        
    protected:
        

        
        /*!
         *    this object, if non-NULL is user-provided to perform actions
         *    after assembly and before returning to the solver
         */
        MAST::WarpingAssembly::PostAssemblyOperation* _post_assembly;

        /*!
         *   L2 norm of the last-assembled residual
         */
        Real _res_l2_norm, _first_iter_res_l2_norm;
        
        /**
         * Defines if Jacobian should be forced to be symmetric. 
         * Default is false.
         */
        bool _force_jacobian_symmetry = false;
        
    };
}


#endif //__mast__warping_assembly__
