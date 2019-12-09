#ifndef __mast__displacement_output__
#define __mast__displacement_output__

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
     *  Retrieves weighted sum of displacements at user-specified DOFs.
     *  This is useful for weighted-sum multiobjected optimization method.
     *  Special case of this is retieving displacement at a single DOF.
     */
    class DisplacementOutput:
    public MAST::OutputAssemblyElemOperations {
        
    public:
        
        /*!
         *    default constructor
         */
        DisplacementOutput(const std::vector<Real> w);
        
        virtual ~DisplacementOutput();
        
        
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
            
            libmesh_error_msg("evaluate_shape_sensitivity not implemented "
                "for displacement output."); // to be implemented
        }
        
        /*!
         * 
         */
        virtual void
        evaluate_topology_sensitivity(const MAST::FunctionBase& f,
                                      const MAST::FieldFunction<RealVectorX>& vel){
            libmesh_error_msg("evalulate_topology_sensitivity not implemented "
                "for displacement output");
        }
        
        /*!
         *   should not get called for this output. Use output_total() instead.
         */
        virtual Real output_for_elem() {
            //
            libmesh_error_msg("output_for_elem() should not get called for "
                "DisplacementOutput. Use output_total() instead.");
        }
        
        /*!
         *   @returns the output quantity value accumulated over all elements
         */
        virtual Real output_total();
        
        /*!
         * 
         */
        virtual Real output_sensitivity_for_elem(const MAST::FunctionBase& p) {
            libmesh_error_msg("output_sensitivity_for_elem not implemented "
                "for DispalcementOutput"); // not yet implemented
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
        
        virtual void evaluate_for_node(const RealVectorX& Xnode,
                                       const RealVectorX& Fpnode) override;
        
        virtual void output_derivative_for_node(const RealVectorX& Xnode, 
            const RealVectorX& Fpnode, RealVectorX& dq_dX) override;
            
        virtual void evaluate_sensitivity_for_node(const MAST::FunctionBase& f,
            const RealVectorX& Xnode, const RealVectorX& dXnode_dparam,
            const RealVectorX& Fpnode, const RealVectorX& dpF_fpparam_node) override;
            
        virtual void get_weights(std::vector<Real>& w);
                
        
    protected:
    
        Real _displacement;
        Real _ddisplacement_dp;
        const std::vector<Real> _w; // Weights for weighted sum of displacements
    };
}

#endif // __mast__displacement_output__
