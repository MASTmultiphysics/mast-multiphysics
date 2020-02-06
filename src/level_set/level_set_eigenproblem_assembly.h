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

#ifndef __mast__level_set_eigenproblem_assembly_h__
#define __mast__level_set_eigenproblem_assembly_h__

// MAST includes
#include "base/eigenproblem_assembly.h"


namespace MAST {
    
    // Forward declerations
    template <typename ValType> class FieldFunction;
    class LevelSetIntersection;
    
    
    class LevelSetEigenproblemAssembly:
    public MAST::EigenproblemAssembly {
    public:
        
        
        /*!
         *   constructor associates this assembly object with the system
         */
        LevelSetEigenproblemAssembly();
        
        
        /*!
         *   destructor resets the association of this assembly object with
         *   the system
         */
        virtual ~LevelSetEigenproblemAssembly();
        
        
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
         *    assembles the matrices for eigenproblem depending on the analysis type
         */
        virtual void
        eigenproblem_assemble(libMesh::SparseMatrix<Real>* A,
                              libMesh::SparseMatrix<Real>* B);
        
        /**
         * Assembly function.  This function will be called
         * to assemble the sensitivity of eigenproblem matrices.
         * The method provides dA/dp and dB/dp for \p f parameter.
         *
         * If the routine is not able to provide sensitivity for this parameter,
         * then it should return false, and the system will attempt to use
         * finite differencing.
         */
        virtual bool
        eigenproblem_sensitivity_assemble (const MAST::FunctionBase& f,
                                           libMesh::SparseMatrix<Real>* sensitivity_A,
                                           libMesh::SparseMatrix<Real>* sensitivity_B);
        
    protected:
        
        MAST::FieldFunction<Real>            *_level_set;
        
        MAST::LevelSetIntersection           *_intersection;
        
        MAST::FieldFunction<RealVectorX>     *_velocity;
    };
}


#endif //__mast__level_set_eigenproblem_assembly_h__


