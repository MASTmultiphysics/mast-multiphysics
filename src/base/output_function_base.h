/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2015  Manav Bhatia
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

#ifndef __mast__output_function_base__
#define __mast__output_function_base__

// C++ includes
#include <vector>

// MAST includes

// libMesh includes
#include "libmesh/point.h"



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
