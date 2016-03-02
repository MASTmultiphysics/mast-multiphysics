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

#ifndef __mast__complex_solver_base_h__
#define __mast__complex_solver_base_h__

// MAST includes
#include "base/mast_data_types.h"

// libMesh includes
#include "libmesh/numeric_vector.h"


namespace MAST {
    
    
    // Forward declerations
    class ComplexAssemblyBase;
    class ElementBase;
    
    
    /*!
     *   uses a Gauss-Siedel method to solve the complex system of equations
     *   for a system.
     */
    class ComplexSolverBase {
        
    public:
        
        /*!
         *  default constructor
         */
        ComplexSolverBase();
        
        
        /*!
         *  destructor
         */
        virtual ~ComplexSolverBase();
        
        
        /*!
         *  sets the assembly object
         */
        void set_assembly(MAST::ComplexAssemblyBase& assembly);
        

        /*!
         *  clears the assembly
         */
        void clear_assembly();
        
        
        /*!
         *  solves the complex system of equations
         */
        virtual void solve();

        
        /*!
         *  solves the complex system of equations using PCFieldSplit
         */
        virtual void solve_pc_fieldsplit();

        
        /*!
         *  solves the complex system of equations using block matrices
         */
        virtual void solve_block_matrix();

        
        /*!
         *  @returns a reference to the real part of the solution
         */
        libMesh::NumericVector<Real>& real_solution();
        
        
        /*!
         *  @returns a constant reference to the real part of the solution
         */
        const libMesh::NumericVector<Real>& real_solution() const;
        
        
        /*!
         *  @returns a reference to the imaginary part of the solution
         */
        libMesh::NumericVector<Real>& imag_solution();

        
        /*!
         *  @returns a constant reference to the imaginary part of the solution
         */
        const libMesh::NumericVector<Real>& imag_solution() const;


        Real tol;
        
        unsigned int max_iters;
        
    protected:
        
        
        /*!
         *   Associated ComplexAssembly object that provides the
         *   element level quantities
         */
        MAST::ComplexAssemblyBase* _assembly;
        
    };
}




#endif // __mast__complex_solver_base_h__

