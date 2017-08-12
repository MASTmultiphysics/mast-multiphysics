/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2017  Manav Bhatia
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

#ifndef __mast__transient_assembly__
#define __mast__transient_assembly__

// MAST includes
#include "base/nonlinear_implicit_assembly.h"



namespace MAST {
    
    // Forward declerations
    class TransientSolverBase;
    
    
    class TransientAssembly:
    public MAST::NonlinearImplicitAssembly {
    public:
        
        /*!
         *   constructor associates this assembly object with the system
         */
        TransientAssembly();
        
        
        /*!
         *   destructor resets the association of this assembly object with
         *   the system
         */
        virtual ~TransientAssembly();
        
        /*!
         *   attaches a system to this discipline, and vice-a-versa
         */
        virtual void
        attach_discipline_and_system(MAST::PhysicsDisciplineBase& discipline,
                                     MAST::SystemInitialization& system) {
            
            libmesh_error_msg("Error! Invalid function call for TransientSystemAssembly.");
        }
        
        /*!
         *   attaches a system to this discipline, and vice-a-versa
         */
        virtual void
        attach_discipline_and_system(MAST::PhysicsDisciplineBase& discipline,
                                     MAST::TransientSolverBase& solver,
                                     MAST::SystemInitialization& sys);
        
        
        /*!
         *   Reattaches to the same system that was attached earlier.
         *
         *   This cannot be called if the clear_discipline_and_system() method
         *   has been called.
         */
        virtual void
        reattach_to_system();
        
        
        /*!
         *   clears association with a system to this discipline, and vice-a-versa
         */
        virtual void
        clear_discipline_and_system();
        
        
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
        
        /**
         * Assembly function.  This function will be called
         * to assemble the sensitivity of system residual prior to a solve and must
         * be provided by the user in a derived class. The method provides dR/dp_i
         * for \par i ^th parameter in the vector \par parameters.
         *
         * If the routine is not able to provide sensitivity for this parameter,
         * then it should return false, and the system will attempt to use
         * finite differencing.
         */
        virtual bool
        sensitivity_assemble (const libMesh::ParameterVector& parameters,
                              const unsigned int i,
                              libMesh::NumericVector<Real>& sensitivity_rhs);
        
        
        
        //**************************************************************
        //these methods are provided for use by the solvers
        //**************************************************************
        
        
        /*!
         *   performs the element calculations over \par elem, for a
         *   system of the form
         *   \f[ f_m(x,\dot{x}) + f_x(x) = 0 \f].
         *   \param f_x     = \f$ f(x) \f$
         *   \param f_m    =  \f$ f_m(x,\dot{x}) \f$
         *   \param f_m_jac_xdot = \f$ \frac{\partial (f_m(x)}{\partial \dot{x}} \f$
         *   \param f_m_jac = \f$ \frac{\partial (f_m(x)}{\partial x} \f$
         *   \param f_x_jac = \f$ \frac{\partial (f(x)}{\partial x} \f$
         *   the element vector and matrix quantities in \par mat and
         *   \par vec, respectively. \par if_jac tells the method to also
         *   assemble the Jacobian, in addition to the residual vector.
         */
        virtual void _elem_calculations(MAST::ElementBase& elem,
                                        bool if_jac,
                                        RealVectorX& f_m,
                                        RealVectorX& f_x,
                                        RealMatrixX& f_m_jac_xdot,
                                        RealMatrixX& f_m_jac,
                                        RealMatrixX& f_x_jac) = 0;
        
        
        
        /*!
         *   performs the element calculations over \par elem, for a
         *   system of the form
         *   \f[ f_m(x,\ddot{x}, \dot{x}) + f_x(x, \dot{x}) = 0 \f].
         *   \param f_x     = \f$ f(x,\dot{x})  \f$
         *   \param f_m    =  \f$ f_m(x,\dot{x}) \f$
         *   \param f_m_jac_xddot = \f$ \frac{\partial (f_m(x)}{\partial \ddot{x}} \f$
         *   \param f_m_jac_xdot = \f$ \frac{\partial (f_m(x)}{\partial \dot{x}} \f$
         *   \param f_m_jac = \f$ \frac{\partial (f_m(x)}{\partial x} \f$
         *   \param f_x_jac_xdot = \f$ \frac{\partial (f(x)}{\partial \dot{x}} \f$
         *   \param f_x_jac = \f$ \frac{\partial (f(x)}{\partial x} \f$
         *   the element vector and matrix quantities in \par mat and
         *   \par vec, respectively. \par if_jac tells the method to also
         *   assemble the Jacobian, in addition to the residual vector.
         */
        virtual void _elem_calculations(MAST::ElementBase& elem,
                                        bool if_jac,
                                        RealVectorX& f_m,
                                        RealVectorX& f_x,
                                        RealMatrixX& f_m_jac_xddot,
                                        RealMatrixX& f_m_jac_xdot,
                                        RealMatrixX& f_m_jac,
                                        RealMatrixX& f_x_jac_xdot,
                                        RealMatrixX& f_x_jac) = 0;
        
        /*!
         *   Calculates the product of Jacobian-solution, and Jacobian-velocity
         *   over the element for a system of the form
         *   \f[ f_m(x,\dot{x}) + f_x(x) = 0 \f].
         *   \param f    = \f$ df(x)/dx \cdot dx +
         *    df_m(x,\dot{x})/dx \cdot dx +
         *    df_m(x,\dot{x})/d\dot{x} \cdot d{\dot x} \f$
         */
        virtual void
        _linearized_jacobian_solution_product(MAST::ElementBase& elem,
                                              RealVectorX& f) = 0;
        
        
        /*!
         *   performs the element sensitivity calculations over \par elem,
         *   and returns the element residual sensitivity in \par vec .
         */
        virtual void _elem_sensitivity_calculations(MAST::ElementBase& elem,
                                                    RealVectorX& vec) = 0;
        
        
    protected:
        
        /*!
         *   performs the element calculations over \par elem, and returns
         *   the element vector and matrix quantities in \par mat and
         *   \par vec, respectively. \par if_jac tells the method to also
         *   assemble the Jacobian, in addition to the residual vector.
         */
        virtual void _elem_calculations(MAST::ElementBase& elem,
                                        bool if_jac,
                                        RealVectorX& vec,
                                        RealMatrixX& mat) {
            
            libmesh_error_msg("Error! Invalid function call for TransientSystemAssembly.");
        }
        
        /*!
         *   performs the element sensitivity calculations over \par elem,
         *   and returns the element residual sensitivity in \par vec .
         */
        virtual void _elem_sensitivity_calculations(MAST::ElementBase& elem,
                                                    bool if_jac,
                                                    RealVectorX& vec,
                                                    RealMatrixX& mat) {
            
            libmesh_error_msg("Error! Invalid function call for TransientSystemAssembly.");
        }
        
        /*!
         *   Pointer to a transient solver object that combines the
         *   element transient data appropriately for the nonlinear solver.
         */
        MAST::TransientSolverBase* _transient_solver;
    };
    
    
}

#endif // __mast__transient_assembly__
