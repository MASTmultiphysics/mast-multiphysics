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

#ifndef __mast__surface_integrated_pressure_output__
#define __mast__surface_integrated_pressure_output__

// C++ includes
#include <map>
#include <vector>

// MAST includes
#include "base/mast_data_types.h"
#include "base/output_assembly_elem_operations.h"


namespace MAST {
    
    
    // Forward declerations
    class FunctionBase;
    
    
    /*!
     *   The surface integrated pressure calculation in the fluid element 
     *   will provide a force vector for the integrated load. This data
     *   structure will provide the mechanism to either select one component 
     *   of the force vector as the output quantity of interest, or the 
     *   magnitude of the force vector. The user must provide the
     */
    class SurfaceIntegratedPressureOutput:
    public MAST::OutputAssemblyElemOperations {
        
    public:
        
        enum OutputMode {
            L2_NORM,     //  L2 norm of the lift
            UNIT_VEC     //  dot product of the lift with a user provided unit vector
        };
        
        
        /*!
         *    default constructor
         */
        SurfaceIntegratedPressureOutput(MAST::SurfaceIntegratedPressureOutput::OutputMode o,
                                        const RealVectorX& n_vec = RealVectorX::Zero(3,1));
        
        virtual ~SurfaceIntegratedPressureOutput();
        
        
        /*!
         *   clears the stored data
         */
        void clear();
        
        /*!
         *   sets the mode that will be used for calculation of the output 
         *   functional from the load vector. If \p o is type \p UNIT_VEC, then
         *   the user must provide the second argument \p n_vec as a 3x1 vector.
         *   The vector will be scaled to unit length.
         */
        void set_output_mode(MAST::SurfaceIntegratedPressureOutput::OutputMode o,
                             const RealVectorX* n_vec = nullptr);
        

        /*!
         *   sets the value of the load
         */
        void set_load(const RealVectorX& v) {
            
            libmesh_assert_equal_to(v.size(), 3);
            _load = v;
        }
        

        /*!
         *    @returns the output functional
         */
        Real value() const;

        
        /*!
         *   sets the value of the load sensitivity wrt function \p f.
         */
        void set_load_sensitivity(const MAST::FunctionBase& f,
                                  const RealVectorX& v) {
            
            libmesh_assert_equal_to(v.size(), 3);
            libmesh_assert(_load_sensitivity.find(&f) == _load_sensitivity.end());
            
            _load_sensitivity[&f] = v;
        }

        
        /*!
         *    @returns the sensitivity of the output functional with respect 
         *    to the provided parameter
         */
        Real sensitivity(const MAST::FunctionBase& f) const;

        
        /*!
         *   sets the value of the load sensitivity wrt function \p f.
         */
        void set_load_derivative(const RealMatrixX& m) {
            
            libmesh_assert_equal_to(m.rows(), 3);
            
            _dload_dX = m;
        }

        /*!
         *    @returns the derivative of the output functional with respect
         *    to the state vector
         */
        RealVectorX derivative() const;
        
        
    protected:

        /*!
         *   output calculation mode
         */
        MAST::SurfaceIntegratedPressureOutput::OutputMode _mode;
        
        
        /*!
         *   unit vector used for output, if \p _mode = \p UNIT_VEC
         */
        RealVectorX _n_vec;
        
        
        /*!
         *   This is the 3x1 vector of the integrated load.
         */
        RealVectorX  _load;
        
        /*!
         *  map of sensitivity of the stress with respect to a parameter
         */
        std::map<const MAST::FunctionBase*, RealVectorX> _load_sensitivity;
        
        
        /*!
         *   derivative of load wrt state vector. The matrix has dimension Nx3
         *   for an element.
         */
        RealMatrixX   _dload_dX;
    };
}

#endif // __mast__surface_integrated_pressure_output__
