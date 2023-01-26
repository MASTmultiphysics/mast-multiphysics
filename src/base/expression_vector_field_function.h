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

#ifndef __mast__expression_vector_field_function__
#define __mast__expression_vector_field_function__

// MAST includes
#include "base/field_function_base.h"

// C++ Mathematical Expression Toolkit Library 
#include "exprtk.h"


namespace MAST {
    
    // Forward declerations
    class Parameter;
    
    class ExpressionVectorFieldFunction:
    public MAST::FieldFunction<RealVectorX> {
    
    public:

        ExpressionVectorFieldFunction(const std::string& nm,
                                      const std::vector<std::string> expr,
                                      const std::set<MAST::Parameter*> params);

        ExpressionVectorFieldFunction(const std::string& nm,
                                      const std::vector<std::string> expr,
                                      const std::map<const std::string, const std::vector<std::string>> dexpr,
                                      const std::set<MAST::Parameter*> params);

        virtual ~ExpressionVectorFieldFunction();

        
        /*!
         *    calculates the value of the function at the specified point,
         *    \p p, and time, \p t, and returns it in \p v.
         */
        virtual void operator() (const libMesh::Point& p,
                                 const Real t,
                                 RealVectorX& v);
        
        
        /*!
         *    calculates the derivative of the function at the specified point,
         *    \p p, and time, \p t, with respect to \p f and returns it in \p v.
         */
        virtual void derivative (const MAST::FunctionBase& f,
                                 const libMesh::Point& p,
                                 const Real t,
                                 RealVectorX& v);


        void derivative (const MAST::FunctionBase& f,
                         const libMesh::Point& p,
                         const Real t,
                         RealVectorX& v,
                         const Real h);

    
        void derivative (const std::string& f,
                         const libMesh::Point& p,
                         const Real t,
                         RealVectorX& v);


        void derivative (const std::string& f,
                         const libMesh::Point& p,
                         const Real t,
                         RealVectorX& v,
                         const Real h);


        const std::vector<std::string> get_expression() const;

        
    protected:

        const std::vector<std::string> _expr;
        const std::map<const std::string, const std::vector<std::string>> _dexpr;
        std::vector<exprtk::expression<Real>> _expression;
        std::map<const std::string, std::vector<exprtk::expression<Real>>> _dexpression;
        Real _x, _y, _z, _t;  // Spatial coordinates (x,y,z) and time (t)
        std::unordered_map<std::string, Real> _param_vars;
        std::unordered_map<std::string, MAST::Parameter*> _params;
    };
}


#endif // __mast__expression_vector_field_function__
