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
         *   initializes the object for the geometric element \p elem. This
         *   expects the object to be in a cleared state, so the user should
         *   call \p clear_elem() between successive initializations.
         */
        virtual void
        init(const libMesh::Elem& elem);
        
        /*!
         *   if set to \p true, the sensitivity calculation will include the
         *   sensitivity of base solution.
         */
        void use_base_sol_for_sensitivity(bool f) {
            _base_sol = f;
        }
        
        /*!
         *   some simulations frequently deal with 1D/2D elements in 3D space,
         *   which requires use of MAST::LocalElemFE.
         */
        virtual bool
        if_use_local_elem() const {
            return true;
        }
        

        /*!
         *   sets additional data for local elem FE.
         */
        virtual void set_local_fe_data(MAST::LocalElemFE& fe,
                                       const libMesh::Elem& e) const;

        virtual void elem_calculations(bool if_jac,
                                       RealVectorX& vec,
                                       RealMatrixX& mat);
        virtual void elem_aerodynamic_force_calculations(ComplexVectorX& vec);
        virtual void elem_sensitivity_calculations(bool if_jac,
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
