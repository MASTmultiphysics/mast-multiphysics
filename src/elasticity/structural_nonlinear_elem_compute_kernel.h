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

#ifndef __mast__structural_nonlinear_assembly_elem_operations__
#define __mast__structural_nonlinear_assembly_elem_operations__

// MAST includes
#include "base/nonlinear_implicit_assembly_elem_operations.h"
#include "elasticity/elasticity_elem_operations.h"

namespace MAST {
    

    // Forward declerations
    class StructuralAssembly;
    
    
    class StructuralNonlinearAssemblyComputations:
    public MAST::NonlinearImplicitAssemblyElemOperations {
        
    public:
        
        /*!
         *   constructor associates this assembly object with the system
         */
        StructuralNonlinearAssemblyComputations():
        MAST::NonlinearImplicitAssemblyElemOperations() { }
        
        
        /*!
         *   destructor resets the association of this assembly object with
         *   the system
         */
        virtual ~StructuralNonlinearAssemblyComputations() {
            
            delete _elem_ops;
        }
        
        /*!
         *   sets the element solution(s) before calculations
         */
        virtual void set_elem_solution(const RealVectorX& sol) {
            
            libmesh_assert_msg(_elem_ops, "Object not initialized");
            
            _elem_ops->compute(vec, if_jac?&mat:nullptr);
        }
        
        /*!
         *   performs the element calculations over \p elem, and returns
         *   the element vector and matrix quantities in \p mat and
         *   \p vec, respectively. \p if_jac tells the method to also
         *   assemble the Jacobian, in addition to the residual vector.
         */
        virtual void elem_calculations(bool if_jac,
                                       RealVectorX& vec,
                                       RealMatrixX& mat) {
            
            libmesh_assert_msg(_elem_ops, "Object not initialized");
            
            _elem_ops->compute(vec, if_jac?&mat:nullptr);
        }
        
        
        /*!
         *   performs the element calculations over \p elem, and returns
         *   the element vector quantity in \p vec. The vector quantity only
         *   include the \f$ [J] \{dX\} \f$ components, so the inherited classes
         *   must ensure that no component of constant forces (traction/body
         *   forces/etc.) are added to this vector.
         */
        virtual void
        elem_linearized_jacobian_solution_product(RealVectorX& vec)
        { libmesh_assert_msg(false, "not implemented yet.");}
        
        
        /*!
         *   performs the element sensitivity calculations over \p elem,
         *   and returns the element residual sensitivity in \p vec .
         */
        virtual void elem_sensitivity_calculations(const MAST::FunctionBase& f,
                                                   RealVectorX& vec) {

            libmesh_assert_msg(_elem_ops, "Object not initialized");
            
            _elem_ops->compute_sensitivity(f, vec);
        }

        
        /*!
         *   performs the element shape sensitivity calculations over \p elem,
         *   and returns the element residual sensitivity in \p vec .
         */
        virtual void
        elem_shape_sensitivity_calculations(const MAST::FunctionBase& f,
                                            RealVectorX& vec)
        { libmesh_assert_msg(false, "not implemented yet.");}

        
        /*!
         *   performs the element topology sensitivity calculations over \p elem,
         *   and returns the element residual sensitivity in \p vec .
         */
        virtual void
        elem_topology_sensitivity_calculations(const MAST::FunctionBase& f,
                                               RealVectorX& vec)
        { libmesh_assert_msg(false, "not implemented yet.");}



        /*!
         *   performs the element topology sensitivity calculations over \p elem,
         *   and returns the element residual sensitivity in \p vec .
         */
        virtual void
        elem_topology_sensitivity_calculations(const MAST::FunctionBase& f,
                                               const MAST::FieldFunction<RealVectorX>& vel,
                                               RealVectorX& vec)
        { libmesh_assert_msg(false, "not implemented yet.");}


        /*!
         *   calculates \f$ d ([J] \{\Delta X\})/ dX  \f$ over \p elem,
         *   and returns the matrix in \p vec .
         */
        virtual void
        elem_second_derivative_dot_solution_assembly(RealMatrixX& mat)
        { libmesh_assert_msg(false, "not implemented yet.");}


        /*!
         *   sets the structural element y-vector if 1D element is used.
         */
        virtual void
        set_elem_data(unsigned int dim,
                      const libMesh::Elem& ref_elem,
                      MAST::GeomElem& elem) const
        { libmesh_assert_msg(false, "not implemented yet.");}


        /*!
         *   initializes the object for the geometric element \p elem. This
         *   expects the object to be in a cleared state, so the user should
         *   call \p clear_elem() between successive initializations.
         */
        virtual void
        init(const MAST::GeomElem& elem) {
            
            libmesh_assert_msg(!_elem_ops, "Object already initialized");
            
            _elem_ops = MAST::ElasticityElemOperations<MAST::ElasticityTraits<Real, Real, Real, MAST::ElasticityElemOperationsContext>>;
            _elem_ops->init(_system->system(), *_discipline);
        }


        /*!
         *   clears the element initialization
         */
        virtual void clear_elem() {
            
            if (_elem_ops) {
                
                delete _elem_ops;
                _elem_ops = nullptr;
            }
            
            MAST::NonlinearImplicitAssemblyElemOperations::clear_elem();
        }
        

    protected:
        
        
        MAST::ElasticityElemOperations<MAST::ElasticityTraits<Real, Real, Real, MAST::ElasticityElemOperationsContext>> *_elem_ops;
    };
}


#endif // __mast__structural_nonlinear_assembly_elem_operations__
