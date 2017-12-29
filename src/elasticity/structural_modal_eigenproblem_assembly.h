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

#ifndef __mast__structural_modal_eigenproblem_assembly__
#define __mast__structural_modal_eigenproblem_assembly__

// MAST includes
#include "base/eigenproblem_assembly.h"


namespace MAST {
    
    
    class StructuralModalEigenproblemAssembly:
    public MAST::EigenproblemAssembly {
    public:
        
        /*!
         *   constructor associates the eigen system with this assembly object
         */
        StructuralModalEigenproblemAssembly();
        
        /*!
         *   destructor resets the association with the eigen system
         *   from this assembly object
         */
        virtual ~StructuralModalEigenproblemAssembly();

        
        /*!
         *   assembles the global A and B matrices for the modal 
         *   eigenvalue problem
         */
        virtual void
        eigenproblem_assemble(libMesh::SparseMatrix<Real>* A,
                              libMesh::SparseMatrix<Real>* B);
        
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
        eigenproblem_sensitivity_assemble (const libMesh::ParameterVector& parameters,
                                           const unsigned int i,
                                           libMesh::SparseMatrix<Real>* sensitivity_A,
                                           libMesh::SparseMatrix<Real>* sensitivity_B);
        

    protected:
        
        /*!
         *   @returns a smart-pointer to a newly created element for
         *   calculation of element quantities.
         */
        virtual std::unique_ptr<MAST::ElementBase>
        _build_elem(const libMesh::Elem& elem);

        /*!
         *   performs the element calculations over \par elem, and returns
         *   the element matrices for the eigenproblem
         *   \f$ A x = \lambda B x \f$.
         */
        virtual void
        _elem_calculations(MAST::ElementBase& elem,
                           RealMatrixX& mat_A,
                           RealMatrixX& mat_B);
        
        /*!
         *   performs the element sensitivity calculations over \par elem,
         *   and returns the element matrices for the eigenproblem
         *   \f$ A x = \lambda B x \f$.
         */
        virtual void
        _elem_sensitivity_calculations(MAST::ElementBase& elem,
                                       RealMatrixX& mat_A,
                                       RealMatrixX& mat_B);
        
        /*!
         *   map of local incompatible mode solution per 3D elements
         */
        std::map<const libMesh::Elem*, RealVectorX> _incompatible_sol;
    };
    
}

#endif // __mast__structural_modal_eigenproblem_assembly__
