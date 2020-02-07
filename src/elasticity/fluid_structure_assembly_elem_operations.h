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

#ifndef __mast__fluid_structure_assembly_elem_operations_h__
#define __mast__fluid_structure_assembly_elem_operations_h__

// MAST includes
#include "base/assembly_elem_operation.h"
#include "elasticity/structural_fluid_interaction_assembly.h"

namespace MAST {
    
    class FluidStructureAssemblyElemOperations:
    public MAST::AssemblyElemOperations {
        
    public:
        
        FluidStructureAssemblyElemOperations();
        virtual ~FluidStructureAssemblyElemOperations();
        
        
        void set_qty_to_evaluate(MAST::StructuralQuantityType q) {
            _qty_type = q;
        }
        
        /*!
         *   sets the structural element y-vector if 1D element is used. 
         */
        virtual void
        set_elem_data(unsigned int dim,
                      const libMesh::Elem& ref_elem,
                      MAST::GeomElem& elem) const;

        /*!
         *   initializes the object for the geometric element \p elem. This
         *   expects the object to be in a cleared state, so the user should
         *   call \p clear_elem() between successive initializations.
         */
        virtual void
        init(const MAST::GeomElem& elem);
        
        /*!
         *   if set to \p true, the sensitivity calculation will include the
         *   sensitivity of base solution.
         */
        void use_base_sol_for_sensitivity(bool f) {
            _base_sol = f;
        }
        
        virtual void elem_calculations(bool if_jac,
                                       RealVectorX& vec,
                                       RealMatrixX& mat);
        virtual void elem_aerodynamic_force_calculations(ComplexVectorX& vec);
        virtual void elem_sensitivity_calculations(const MAST::FunctionBase& f,
                                                   bool if_jac,
                                                   RealVectorX& vec,
                                                   RealMatrixX& mat);
        virtual void elem_second_derivative_dot_solution_assembly(RealMatrixX& m);
    protected:
        
        /*!
         *   whether or not the base solution was included for linearization
         */
        bool _base_sol;
        
        /*!
         *   this defines the quantity to be assembled
         */
        MAST::StructuralQuantityType _qty_type;
    };
}



#endif // __mast__fluid_structure_assembly_elem_operations_h__
