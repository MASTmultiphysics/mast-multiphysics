
#ifndef __mast__output_function_base__
#define __mast__output_function_base__

// C++ includes
#include <vector>

// MAST includes

// libMesh includes
#include "geom/point.h"



namespace MAST {
    
    // identifies the evaluation mode for pointwise quantities
    enum PointwiseOutputEvaluationMode {
        CENTROID,    // at the element centroid
        ELEM_QP,     // at the element quadrature points
        SPECIFIED_QP // at quadrature points specified by the user
    };
    
    
    //template <typename ValType>
    class OutputFunctionBase {
        
    public:
        
        OutputFunctionBase();

        /*!
         *    virtual destructor
         */
        virtual ~OutputFunctionBase();
       
        
        /*!
         *   sets the points of evaluation
         */
        void set_qp_for_evaluation(const std::vector<libMesh::Point>& pts);

        
        /*!
         *   returns the points of evaluation
         */
        const std::vector<libMesh::Point>& get_qp_for_evaluation() const;

        
        /*!
         *   @returns the mode of output evaluation
         */
        MAST::PointwiseOutputEvaluationMode evaluation_mode() const {
            
            return _eval_mode;
        }
        

        
        
    protected:
        

        std::vector<libMesh::Point> _eval_points;
        
        MAST::PointwiseOutputEvaluationMode _eval_mode;
        
    };
}

#endif // __mast__output_function_base__
