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

#ifndef __mast__complex_mesh_field_function__
#define __mast__complex_mesh_field_function__

// MAST includes
#include "base/field_function_base.h"


// libMesh includes
#include "libmesh/numeric_vector.h"
#include "libmesh/mesh_function.h"



namespace MAST {
    
    // Forward declerations
    class SystemInitialization;
    
    
    /*!
     *    This provides a wrapper FieldFunction compatible class that
     *    interpolates the solution using libMesh's MeshFunction class.
     */
    class ComplexMeshFieldFunction:
    public MAST::FieldFunction<ComplexVectorX> {
        
    public:
        /*!
         *   constructor
         */
        ComplexMeshFieldFunction(MAST::SystemInitialization& sys,
                                 const std::string& nm);
        
        
        /*!
         *   destructor
         */
        virtual ~ComplexMeshFieldFunction();
        
        
        virtual void operator() (const libMesh::Point& p,
                                 const Real t,
                                 ComplexVectorX& v) const;
        
        
        virtual void perturbation (const libMesh::Point& p,
                                   const Real t,
                                   ComplexVectorX& v) const;
        
        
        void init(const libMesh::NumericVector<Real>& sol_re,
                  const libMesh::NumericVector<Real>& sol_im);
        
        
        void init_perturbation(const libMesh::NumericVector<Real>& dsol_re,
                               const libMesh::NumericVector<Real>& dsol_im);
        
        
        /*!
         *    @returns a reference to the libMesh mesh function
         */
        std::pair<libMesh::MeshFunction*, libMesh::MeshFunction*>
        get_function() {
            
            libmesh_assert(_function_re);
            
            return std::pair<libMesh::MeshFunction*, libMesh::MeshFunction*>
            (_function_re, _function_im);
        }
        
        /*!
         *    @returns a reference to the libMesh mesh function for the
         *    perturbation in solution
         */
        std::pair<libMesh::MeshFunction*, libMesh::MeshFunction*>
        get_perturbed_function() {
            
            libmesh_assert(_perturbed_function_re);
            
            return std::pair<libMesh::MeshFunction*, libMesh::MeshFunction*>
            (_perturbed_function_re, _perturbed_function_im);
        }
        
        
        /*!
         *   clear the solution and mesh function data structures
         */
        void clear();
        
    protected:
        
        /*!
         *  current system for which solution is to be interpolated
         */
        MAST::SystemInitialization* _system;
        
        /*!
         *   current solution that is going to be interpolated
         */
        libMesh::NumericVector<Real>
        *_sol_re,
        *_sol_im,
        *_perturbed_sol_re,
        *_perturbed_sol_im;
        
        /*!
         *   the MeshFunction object that performs the interpolation
         */
        libMesh::MeshFunction
        *_function_re,
        *_function_im,
        *_perturbed_function_re,
        *_perturbed_function_im;
        
    };
}

#endif // __mast__complex_mesh_field_function__

