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
    
    // identifies the various output quantities implemented in the code
    enum OutputQuantityType {
        
        STRAIN_STRESS_TENSOR
    };
    
    
    // identifies the evaluation mode for pointwise quantities
    enum PointwiseOutputEvaluationMode {
        CENTROID,    // at the element centroid
        ELEM_QP,     // at the element quadrature points
        SPECIFIED_POINTS // at quadrature points specified by the user
    };
    
    
    //template <typename ValType>
    class OutputFunctionBase {
        
    public:
        
        OutputFunctionBase(MAST::OutputQuantityType t);

        /*!
         *    virtual destructor
         */
        virtual ~OutputFunctionBase();
       
       
        /*!
         *   @returns the type of this output quantity
         */
        MAST::OutputQuantityType type() const;
        
        
        /*!
         *   sets the quadrature points of evaluation of the quantity
         */
        void set_points_for_evaluation(const std::vector<libMesh::Point>& pts);

        
        /*!
         *   returns the points of evaluation
         */
        const std::vector<libMesh::Point>& get_points_for_evaluation() const;

        
        /*!
         *   sets the mode of output evaluation
         */
        void set_evaluation_mode(MAST::PointwiseOutputEvaluationMode m)  {
            
            _eval_mode = m;
        }

        
        
        /*!
         *   @returns the mode of output evaluation
         */
        MAST::PointwiseOutputEvaluationMode evaluation_mode() const {
            
            return _eval_mode;
        }
        

        
        
    protected:

        
        /*!
         *   Type of this quantity
         */
        MAST::OutputQuantityType _type;
        

        std::vector<libMesh::Point> _eval_points;
        
        MAST::PointwiseOutputEvaluationMode _eval_mode;
        
    };
}

#endif // __mast__output_function_base__
