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
    class Parameter;
    
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
         *   @returns a reference to the assembly object
         */
        inline
        MAST::ComplexAssemblyBase& get_assembly() {
            
            return *_assembly;
        }

        
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
         *  solves the complex system of equations using block matrices. If 
         *  no argument is specified for \par p, then the system is solved. 
         *  Otherwise, the sensitivity of the system is solved with respect
         *  to the parameter p
         */
        virtual void solve_block_matrix(MAST::Parameter* p = nullptr);

        
        /*!
         *  @returns a reference to the real part of the solution. If 
         *  \par if_sens is true, the the sensitivity vector is returned. Note,
         *  that the sensitivity can be requested only after a sensitivity 
         *  solve.
         */
        libMesh::NumericVector<Real>& real_solution(bool if_sens=false);
        
        
        /*!
         *  @returns a constant reference to the real part of the solution. If
         *  \par if_sens is true, the the sensitivity vector is returned. Note,
         *  that the sensitivity can be requested only after a sensitivity
         *  solve.
         */
        const libMesh::NumericVector<Real>& real_solution(bool if_sens=false) const;
        
        
        /*!
         *  @returns a reference to the imaginary part of the solution. If
         *  \par if_sens is true, the the sensitivity vector is returned. Note,
         *  that the sensitivity can be requested only after a sensitivity
         *  solve.
         */
        libMesh::NumericVector<Real>& imag_solution(bool if_sens=false);

        
        /*!
         *  @returns a constant reference to the imaginary part of the solution.
         *  If \par if_sens is true, the the sensitivity vector is returned.
         *  Note, that the sensitivity can be requested only after a sensitivity
         *  solve.
         */
        const libMesh::NumericVector<Real>& imag_solution(bool if_sens=false) const;


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

