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

#ifndef __mast__level_set_nonlinear_implicit_assembly_h__
#define __mast__level_set_nonlinear_implicit_assembly_h__

// MAST includes
#include "base/nonlinear_implicit_assembly.h"


namespace MAST {
    
    // Forward declerations
    template <typename ValType> class FieldFunction;
    class LevelSetIntersection;
    class LevelSetInterfaceDofHandler;
    class LevelSetVoidSolution;
    class FilterBase;
    
    
    class LevelSetNonlinearImplicitAssembly:
    public MAST::NonlinearImplicitAssembly {
    public:
        
        
        /*!
         *   constructor associates this assembly object with the system
         */
        LevelSetNonlinearImplicitAssembly(bool enable_dof_handler);
        
        
        /*!
         *   destructor resets the association of this assembly object with
         *   the system
         */
        virtual ~LevelSetNonlinearImplicitAssembly();

        /*!
         *  sets the flag on whether or not to evaluate the output on
         *  negative level set function
         */
        void set_evaluate_output_on_negative_phi(bool f);

        
        /*!
         *   @return flag if using dof_handler or not
         */
        bool if_use_dof_handler() const;
        
        
        /*!
         *   attaches level set function to \p this
         */
        virtual void
        set_level_set_function(MAST::FieldFunction<Real>& level_set,
                               const MAST::FilterBase& filter);

        
        /*!
         *   attaches indicator function to \p this.
         */
        virtual void
        set_indicator_function(MAST::FieldFunction<RealVectorX>& indicator);
        
        /*!
         *   clears association with level set function
         */
        virtual void
        clear_level_set_function();

        /*!
         *   the velocity function used to calculate topology sensitivity
         */
        virtual void
        set_level_set_velocity_function(MAST::FieldFunction<RealVectorX>& velocity);
        
        
        /*!
         *   clears the velocity function
         */
        virtual void
        clear_level_set_velocity_function();

        
        /*!
         *  @returns a reference to the level set function
         */
        MAST::LevelSetIntersection& get_intersection();

        
        /*!
         *  @returns a reference to the \p LevelSetInterfaceDofHandler object
         */
        MAST::LevelSetInterfaceDofHandler& get_dof_handler();

        
        /*!
         *    function that assembles the matrices and vectors quantities for
         *    nonlinear solution
         */
        virtual void
        residual_and_jacobian (const libMesh::NumericVector<Real>& X,
                               libMesh::NumericVector<Real>* R,
                               libMesh::SparseMatrix<Real>*  J,
                               libMesh::NonlinearImplicitSystem& S);

        
        virtual bool
        sensitivity_assemble (const libMesh::NumericVector<Real>& X,
                              bool if_localize_sol,
                              const MAST::FunctionBase& f,
                              libMesh::NumericVector<Real>& sensitivity_rhs,
                              bool close_vector = true);
        
        
        virtual void
        calculate_output_derivative(const libMesh::NumericVector<Real>& X,
                                    bool if_localize_sol,
                                    MAST::OutputAssemblyElemOperations& output,
                                    libMesh::NumericVector<Real>& dq_dX);
        
//#define MAST_ENABLE_PLPLOT 1
#if MAST_ENABLE_PLPLOT == 1
        void plot_sub_elems(bool plot_reference_elem,
                            bool plot_low_phi_elem,
                            bool plot_high_phi_elem);
#endif
        
        
        /*!
         *   calculates the value of quantity \f$ q(X,p) \f$.
         */
        virtual void
        calculate_output(const libMesh::NumericVector<Real>& X,
                         bool if_localize_sol,
                         MAST::OutputAssemblyElemOperations& output);
        
        
        
        /*!
         *   evaluates the sensitivity of the outputs in the attached
         *   discipline with respect to the parametrs in \p params.
         *   The base solution should be provided in \p X. If total sensitivity
         *   is desired, then \p dXdp should contain the sensitivity of
         *   solution wrt the parameter \p p, otherwise it can be set to
         *   nullptr. If \p dXdp is nullptr,
         *   the calculated sensitivity will be the partial derivarive of
         *   \p output wrt \p p.
         */
        virtual void
        calculate_output_direct_sensitivity(const libMesh::NumericVector<Real>& X,
                                            bool if_localize_sol,
                                            const libMesh::NumericVector<Real>* dXdp,
                                            bool if_localize_sol_sens,
                                            const MAST::FunctionBase& p,
                                            MAST::OutputAssemblyElemOperations& output);
        
        
    protected:

        bool                                  _enable_dof_handler;
        
        bool                                  _evaluate_output_on_negative_phi;

        MAST::FieldFunction<Real>            *_level_set;

        MAST::FieldFunction<RealVectorX>     *_indicator;

        MAST::LevelSetIntersection           *_intersection;

        MAST::LevelSetInterfaceDofHandler    *_dof_handler;
        
        MAST::LevelSetVoidSolution           *_void_solution_monitor;

        MAST::FieldFunction<RealVectorX>     *_velocity;
        
        const MAST::FilterBase               *_filter;

    };
}


#endif //__mast__level_set_nonlinear_implicit_assembly_h__

