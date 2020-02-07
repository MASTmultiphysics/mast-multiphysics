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

#ifndef __mast__compliance_output__
#define __mast__compliance_output__

// C++ includes
#include <map>
#include <vector>

// MAST includes
#include "base/mast_data_types.h"
#include "base/physics_discipline_base.h"
#include "base/output_assembly_elem_operations.h"


// libMesh includes
#include "libmesh/elem.h"

namespace MAST {
    
    
    // Forward declerations
    class FunctionBase;
    
    
    /*!
     *   Computes the compliance as \f$ \frac{1}{2} U^T K U \f$.
     */
    class ComplianceOutput:
    public MAST::OutputAssemblyElemOperations {
        
    public:
        
        /*!
         *    default constructor
         */
        ComplianceOutput();
        
        virtual ~ComplianceOutput();
        
        
        /*!
         *   sets the structural element y-vector if 1D element is used.
         */
        virtual void
        set_elem_data(unsigned int dim,
                      const libMesh::Elem& ref_elem,
                      MAST::GeomElem& elem) const;
        
        /*!
         *   initialize for the element.
         */
        virtual void init(const MAST::GeomElem& elem);
        
        /*!
         *   zeroes the output quantity values stored inside this object
         *   so that assembly process can begin. This will zero out data
         *   so that it is ready for a new evaluation. Before sensitivity
         *   analysis, call the other method, since some nonlinear
         *   functionals need the forward quantities for sensitivity analysis,
         *   eg., stress output.
         */
        virtual void zero_for_analysis();
        
        
        /*!
         *   zeroes the output quantity values stored inside this object
         *   so that assembly process can begin. This will only zero the
         *   data to compute new sensitivity analysis.
         */
        virtual void zero_for_sensitivity();
        
        /*!
         *    this evaluates all relevant stress components on the element to
         *    evaluate the p-averaged quantity.
         *    This is only done on the current element for which this
         *    object has been initialized.
         */
        virtual void evaluate();
        
        /*!
         *    this evaluates all relevant stress sensitivity components on
         *    the element to evaluate the p-averaged quantity sensitivity.
         *    This is only done on the current element for which this
         *    object has been initialized.
         */
        virtual void evaluate_sensitivity(const MAST::FunctionBase& f);
        
        /*!
         *    this evaluates all relevant shape sensitivity components on
         *    the element.
         *    This is only done on the current element for which this
         *    object has been initialized.
         */
        virtual void evaluate_shape_sensitivity(const MAST::FunctionBase& f) {
            
            libmesh_assert(false); // to be implemented
        }
        
        /*!
         *    this evaluates all relevant topological sensitivity components on
         *    the element.
         *    This is only done on the current element for which this
         *    object has been initialized.
         */
        virtual void
        evaluate_topology_sensitivity(const MAST::FunctionBase& f);

        /*!
         *    This evaluates the contribution to the topology sensitivity on the
         *    boundary. Given that the integral is nonlinear due to the \f$p-\f$norm,
         *    the expression is quite involved:
         *   \f[  \frac{ \frac{1}{p} \left( \int_\Omega (\sigma_{VM}(\Omega))^p ~
         *    d\Omega \right)^{\frac{1}{p}-1}}{\left(  \int_\Omega ~ d\Omega \right)^{\frac{1}{p}}}
         *    \int_\Gamma  V_n \sigma_{VM}^p  ~d\Gamma +
         *    \frac{ \frac{-1}{p} \left( \int_\Omega (\sigma_{VM}(\Omega))^p ~
         *    d\Omega \right)^{\frac{1}{p}}}{\left(  \int_\Omega ~ d\Omega \right)^{\frac{1+p}{p}}}
         *    \int_\Gamma  V_n ~d\Gamma \f]
         */
        virtual void
        evaluate_topology_sensitivity(const MAST::FunctionBase& f,
                                      const MAST::FieldFunction<RealVectorX>& vel);
        
        /*!
         *   should not get called for this output. Use output_total() instead.
         */
        virtual Real output_for_elem() {
            //
            libmesh_error();
        }
        
        /*!
         *   @returns the output quantity value accumulated over all elements
         */
        virtual Real output_total();
        
        /*!
         *   @returns the sensitivity of p-norm von Mises stress for the
         *   \f$p-\f$norm identified using \p set_p_val(). The returned quantity
         *   is evaluated for the element for which this object is initialized.
         */
        virtual Real output_sensitivity_for_elem(const MAST::FunctionBase& p) {
            libmesh_error(); // not yet implemented
        }
        
        /*!
         *   @returns the output quantity sensitivity for parameter.
         *   This method calculates the partial derivative of quantity
         *    \f[ \frac{\partial q(X, p)}{\partial p} \f]  with
         *    respect to parameter \f$ p \f$. This returns the quantity
         *   accumulated over all elements.
         */
        virtual Real output_sensitivity_total(const MAST::FunctionBase& p);
        
        
        /*!
         *   calculates the derivative of p-norm von Mises stress for the
         *   \f$p-\f$norm identified using \p set_p_val(). The quantity is
         *   evaluated over the current element for which this object
         *   is initialized.
         */
        virtual void output_derivative_for_elem(RealVectorX& dq_dX);
        
        
    protected:
    
        Real _compliance;
        Real _dcompliance_dp;
    };
}

#endif // __mast__compliance_output__
