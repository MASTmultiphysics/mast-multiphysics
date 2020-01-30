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

#ifndef __mast__mesh_field_function__
#define __mast__mesh_field_function__

// MAST includes
#include "base/field_function_base.h"


// libMesh includes
#include "libmesh/numeric_vector.h"
#include "libmesh/mesh_function.h"
#include "libmesh/system.h"



namespace MAST {

    // Forward declerations
    class SystemInitialization;
    
    
    /*!
     *    This provides a wrapper FieldFunction compatible class that
     *    interpolates the solution using libMesh's MeshFunction class.
     */
    class MeshFieldFunction:
    public MAST::FieldFunction<RealVectorX> {
        
    public:
        /*!
         *   constructor
         */
        MeshFieldFunction(MAST::SystemInitialization& sys,
                          const std::string& nm);

        
        /*!
         *   constructor
         */
        MeshFieldFunction(libMesh::System& sys, const std::string& nm);

        
        /*!
         *   destructor
         */
        virtual ~MeshFieldFunction();

       
        /*!
         *    calculates the value of the function at the specified point,
         *    \p p, and time, \p t, and returns it in \p v.
         */
        virtual void operator() (const libMesh::Point& p,
                                 const Real t,
                                 RealVectorX& v) const;
        
        /*!
         *    calculates the gradient of value of the function at the specified
         *    point, \p p, and time, \p t, and returns it in \p g.
         *    g(i,j) = dv(i)/dx(j)
         */
        virtual void gradient(const libMesh::Point& p,
                              const Real t,
                              RealMatrixX& g) const;

        /*!
         *    calculates the value of perturbation in the function at
         *    the specified point, \p p, and time, \p t, and returns it
         *    in \p v.
         */
        virtual void perturbation (const libMesh::Point& p,
                                   const Real t,
                                   RealVectorX& v) const;

        
        /*!
         *    calculates the value of the function at the specified point,
         *    \p p, and time, \p t, and returns it in \p v.
         */
        virtual void derivative (const MAST::FunctionBase& f,
                                 const libMesh::Point& p,
                                 const Real t,
                                 RealVectorX& v) const;


        /*!
         *    calculates the value of the function at the specified point,
         *    \p p, and time, \p t, and returns it in \p v.
         */
        virtual void perturbation_gradient (const libMesh::Point& p,
                                            const Real t,
                                            RealMatrixX& v) const;

        
        /*!
         *   initializes the data structures to perform the interpolation 
         *   function of \p sol. If \p dsol is provided, then it is used
         *   as the perturbation of \p sol.
         */
        void init(const libMesh::NumericVector<Real>& sol,
                  const libMesh::NumericVector<Real>* dsol = nullptr);

        
        /*!
         *    @returns a reference to the libMesh mesh function
         */
        libMesh::MeshFunction& get_function() {
            
            libmesh_assert(_function);
            return *_function;
        }

        /*!
         *    @returns a reference to the libMesh mesh function for the 
         *    perturbation in solution
         */
        libMesh::MeshFunction& get_perturbed_function() {
            
            libmesh_assert(_perturbed_function);
            return *_perturbed_function;
        }

        
        /*!
         *    When a mesh field function is attached to an assembly routine
         *    during system assembly, then the current solution can be 
         *    provided by the element quadrature point update. This method
         *    allows the element to provide the solution and in turn
         *    override the mesh function evaluation.
         */
        virtual void set_element_quadrature_point_solution(RealVectorX& sol);

        
        /*!
         *    clears the quadrature point solution provided by the 
         *    corresponding set method above.
         */
        virtual void clear_element_quadrature_point_solution();

        
        /*!
         *   clear the solution
         */
        void clear();

    protected:

        /*!
         *  flag is set to true when the quadrature point solution is 
         *  provided by an element
         */
        bool _use_qp_sol;
        
        
        /*!
         *   quadrature point solution of the element
         */
        RealVectorX _qp_sol;
        
        /*!
         *  current system for which solution is to be interpolated
         */
        libMesh::System* _sys;
        
        /*!
         *   current solution that is going to be interpolated
         */
        libMesh::NumericVector<Real> *_sol, *_dsol;
        
        /*!
         *   the MeshFunction object that performs the interpolation
         */
        libMesh::MeshFunction *_function, *_perturbed_function;
    };
}

#endif // __mast__mesh_field_function__

