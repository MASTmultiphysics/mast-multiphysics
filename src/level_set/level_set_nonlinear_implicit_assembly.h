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

#ifndef __mast__level_set_nonlinear_implicit_assembly_h__
#define __mast__level_set_nonlinear_implicit_assembly_h__

// MAST includes
#include "base/nonlinear_implicit_assembly.h"


namespace MAST {
    
    // Forward declerations
    template <typename ValType> class FieldFunction;
    class LevelSetIntersection;
    
    
    class LevelSetNonlinearImplicitAssembly:
    public MAST::NonlinearImplicitAssembly,
    public libMesh::System::Constraint {
    public:
        
        
        /*!
         *   constructor associates this assembly object with the system
         */
        LevelSetNonlinearImplicitAssembly();
        
        
        /*!
         *   destructor resets the association of this assembly object with
         *   the system
         */
        virtual ~LevelSetNonlinearImplicitAssembly();
        
        
        /*!
         *   attaches level set function to \p this
         */
        virtual void
        set_level_set_function(MAST::FieldFunction<Real>& level_set);

        /*!
         *   clears association with level set function
         */
        virtual void
        clear_level_set_function();
        
        /*!
         *  @returns a reference to the level set function
         */
        MAST::LevelSetIntersection& get_intersection();
        
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
         *   calculates the value of quantity \f$ q(X,p) \f$.
         */
        virtual void
        calculate_output(const libMesh::NumericVector<Real>& X,
                         MAST::OutputAssemblyElemOperations& output);
        
        
        
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
         *   provides implementation of the libMesh::System::Constraint::constrain()
         *   virtual method
         */
        virtual void
        constrain ();
        
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

        void _initialize_constraints();
        
        MAST::FieldFunction<Real> *_level_set;

        MAST::LevelSetIntersection *_intersection;

    };
}


#endif //__mast__level_set_nonlinear_implicit_assembly_h__

