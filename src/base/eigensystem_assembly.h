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


#ifndef __mast__eigensystem_assembly_h__
#define __mast__eigensystem_assembly_h__


// libMesh includes
#include "libmesh/sparse_matrix.h"
#include "libmesh/parameter_vector.h"


namespace MAST {
    
    class EigenSystemAssembly {
        
    public:

        /*!
         *   constructor
         */
        EigenSystemAssembly() { }
        
        
        /*!
         *   destructor
         */
        virtual ~EigenSystemAssembly() { }
        
        
        /*!
         *    assembles the matrices for eigenproblem depending on the analysis type
         */
        virtual void
        eigenproblem_assemble(libMesh::SparseMatrix<Real>* matrix_A,
                              libMesh::SparseMatrix<Real>* matrix_B) = 0;
        
        
        /**
         * Assembly function.  This function will be called
         * to assemble the sensitivity of eigenproblem matrices.
         * The method provides dA/dp_i and dB/dpi for \par i ^th parameter
         * in the vector \par parameters.
         *
         * If the routine is not able to provide sensitivity for this parameter,
         * then it should return false, and the system will attempt to use
         * finite differencing.
         */
        virtual bool
        eigenproblem_sensitivity_assemble(const libMesh::ParameterVector& parameters,
                                          const unsigned int i,
                                          libMesh::SparseMatrix<Real>* sensitivity_A,
                                          libMesh::SparseMatrix<Real>* sensitivity_B) = 0;
        
        
    };
}


#endif // __mast__eigensystem_assembly_h__
