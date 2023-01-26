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

// MAST includes
#include "base/expression_vector_field_function.h"
#include "base/parameter.h"


MAST::ExpressionVectorFieldFunction::
ExpressionVectorFieldFunction(const std::string& nm,
                              const std::vector<std::string> expr,
                              const std::set<MAST::Parameter*> params):
MAST::FieldFunction<RealVectorX>(nm), 
_expr(expr){

    exprtk::symbol_table<Real> symbol_table;
    exprtk::parser<Real> parser;

    // Add spatial variables
    symbol_table.add_variable("x", _x);
    symbol_table.add_variable("y", _y);
    symbol_table.add_variable("z", _z);

    // Add time variable
    symbol_table.add_variable("t", _t);

    // Add MAST Parameter variables
    std::set<MAST::Parameter*>::iterator it;
    for (it = params.begin(); it != params.end(); ++it) {
        std::string name = (**it).name();
        _param_vars[name] = (**it)();
        _params[name] = *it;
        symbol_table.add_variable(name, _param_vars[name]);
    }

    uint n_components = _expr.size();

    _expression.resize(n_components);

    for (uint i=0; i<n_components; i++) {
        _expression[i].register_symbol_table(symbol_table);
        if (!(parser.compile(_expr[i], _expression[i])))
        {
            throw std::runtime_error("Error encoutered in exprtk when trying to compile expression.");
        }
    }
}


MAST::ExpressionVectorFieldFunction::
ExpressionVectorFieldFunction(const std::string& nm,
                              const std::vector<std::string> expr,
                              const std::map<const std::string, const std::vector<std::string>> dexpr,
                              const std::set<MAST::Parameter*> params):
MAST::FieldFunction<RealVectorX>(nm), 
_expr(expr),
_dexpr(dexpr){

    exprtk::symbol_table<Real> symbol_table;
    exprtk::parser<Real> parser;

    // Add spatial variables
    symbol_table.add_variable("x", _x);
    symbol_table.add_variable("y", _y);
    symbol_table.add_variable("z", _z);

    // Add time variable
    symbol_table.add_variable("t", _t);

    // Add MAST Parameter variables
    std::set<MAST::Parameter*>::iterator it;
    for (it = params.begin(); it != params.end(); ++it) {
        std::string name = (**it).name();
        _param_vars[name] = (**it)();
        _params[name] = *it;
        symbol_table.add_variable(name, _param_vars[name]);
    }

    uint n_components = _expr.size();
    _expression.resize(n_components);

    for (uint i=0; i<n_components; i++) {
        _expression[i].register_symbol_table(symbol_table);
        if (!(parser.compile(_expr[i], _expression[i])))
        {
            throw std::runtime_error("Error encoutered in exprtk when trying to compile expression.");
        }
    }

    std::map<const std::string, const std::vector<std::string>>::const_iterator it_expr;
    for (it_expr = _dexpr.begin(); it_expr != _dexpr.end(); it_expr++) {

        const std::string var = it_expr->first;
        const std::vector<std::string> dexpr_dvar = it_expr->second;

        uint n_components = _expr.size();
        exprtk::expression<Real> dexpr;
        std::vector<exprtk::expression<Real>> dexpression(n_components, dexpr);
        _dexpression.insert({var, dexpression});

        for (uint i=0; i<n_components; i++) {
            _dexpression.at(var)[i].register_symbol_table(symbol_table);
            if (!(parser.compile(dexpr_dvar[i], _dexpression.at(var)[i])))
            {
                throw std::runtime_error("Error encoutered in exprtk when trying to compile dexpression.");
            }
        }
    }
}


MAST::ExpressionVectorFieldFunction::~ExpressionVectorFieldFunction() {
    
}


void
MAST::ExpressionVectorFieldFunction::operator() (const libMesh::Point& p, 
                                           const Real t, 
                                           RealVectorX& v) {
    // Evaluate spatial coordinates
    _x = p(0);
    _y = p(1);
    _z = p(2);
    
    // Evaluate time
    _t = t;

    // Evaluate Parameters (in case they changed since the constructor)
    for (auto it=_params.begin(); it!=_params.end(); it++) {
        _param_vars[it->first] = (*(it->second))();
    }

    uint n_components = _expression.size();
    for (uint i=0; i<n_components; i++) {
        v[i] = _expression[i].value();
    }
}


void
MAST::ExpressionVectorFieldFunction::derivative (const MAST::FunctionBase& f,
                                           const libMesh::Point& p,
                                           const Real t,
                                           RealVectorX& v) {
    // Evaluate spatial coordinates
    _x = p(0);
    _y = p(1);
    _z = p(2);
    
    // Evaluate time
    _t = t;

    // Evaluate Parameters (in case they changed since the constructor)
    for (auto it=_params.begin(); it!=_params.end(); it++) {
        _param_vars[it->first] = (*(it->second))();
    }

    uint n_components = _expression.size();
    // The symbol table is the same for all expressions so just grab it from any one of them
    exprtk::symbol_table<Real>& symbol_table = _dexpression.begin()->second[0].get_symbol_table();

    if (symbol_table.is_variable(f.name())) {
        // If the parameter we are taking the derivative with respect to (f) 
        // is a varaible, calculate the derivative.
        if (_dexpression.count(f.name())) {
            // Use the user-supplied analytical derivative if available
            for (uint i=0; i<n_components; i++) {
                v[i] = _dexpression.at(f.name())[i].value();
            }
        }
        else {
            // Use exprtk's finite difference approximation if analytical
            // derivative was not available.
            for (uint i=0; i<n_components; i++) {
                v[i] = exprtk::derivative(_expression[i], f.name());
            }
        }
    }
    else {
        // If the parameter we are taking the derivative with respect to (f)
        // is not a varaible, then the derivative is zero.
        for (uint i=0; i<n_components; i++) {
            v[i] = 0.0;
        }
    }
}


void
MAST::ExpressionVectorFieldFunction::derivative (const std::string& f,
                                           const libMesh::Point& p,
                                           const Real t,
                                           RealVectorX& v) {
    // Evaluate spatial coordinates
    _x = p(0);
    _y = p(1);
    _z = p(2);
    
    // Evaluate time
    _t = t;

    // Evaluate Parameters (in case they changed since the constructor)
    for (auto it=_params.begin(); it!=_params.end(); it++) {
        _param_vars[it->first] = (*(it->second))();
    }

    uint n_components = _expression.size();
    // The symbol table is the same for all expressions so just grab it from any one of them
    exprtk::symbol_table<Real>& symbol_table = _dexpression.begin()->second[0].get_symbol_table();

    if (symbol_table.is_variable(f)) {
        // If the parameter we are taking the derivative with respect to (f) 
        // is a varaible, calculate the derivative.
        if (_dexpression.count(f)) {
            // Use the user-supplied analytical derivative if available
            for (uint i=0; i<n_components; i++) {
                v[i] = _dexpression.at(f)[i].value();
            }
        }
        else {
            // Use exprtk's finite difference approximation if analytical
            // derivative was not available.
            for (uint i=0; i<n_components; i++) {
                v[i] = exprtk::derivative(_expression[i], f);
            }
        }
    }
    else {
        // If the parameter we are taking the derivative with respect to (f)
        // is not a varaible, then the derivative is zero.
        for (uint i=0; i<n_components; i++) {
            v[i] = 0.0;
        }
    }
}


void
MAST::ExpressionVectorFieldFunction::derivative (const MAST::FunctionBase& f,
                                           const libMesh::Point& p,
                                           const Real t,
                                           RealVectorX& v,
                                           const Real h) {
    // Evaluate spatial coordinates
    _x = p(0);
    _y = p(1);
    _z = p(2);
    
    // Evaluate time
    _t = t;

    // Evaluate Parameters (in case they changed since the constructor)
    for (auto it=_params.begin(); it!=_params.end(); it++) {
        _param_vars[it->first] = (*(it->second))();
    }

    uint n_components = _expression.size();
    // The symbol table is the same for all expressions so just grab it from any one of them
    exprtk::symbol_table<Real>& symbol_table = _dexpression.begin()->second[0].get_symbol_table();

    if (symbol_table.is_variable(f.name())) {
        // If the parameter we are taking the derivative with respect to (f) 
        // is a varaible, calculate the derivative.
        if (_dexpression.count(f.name())) {
            // Use the user-supplied analytical derivative if available
            for (uint i=0; i<n_components; i++) {
                v[i] = _dexpression.at(f.name())[i].value();
            }
        }
        else {
            // Use exprtk's finite difference approximation if analytical
            // derivative was not available.
            for (uint i=0; i<n_components; i++) {
                v[i] = exprtk::derivative(_expression[i], f.name(), h);
            }
        }
    }
    else {
        // If the parameter we are taking the derivative with respect to (f)
        // is not a varaible, then the derivative is zero.
        for (uint i=0; i<n_components; i++) {
            v[i] = 0.0;
        }
    }
}


void
MAST::ExpressionVectorFieldFunction::derivative (const std::string& f,
                                           const libMesh::Point& p,
                                           const Real t,
                                           RealVectorX& v,
                                           const Real h) {
    // Evaluate spatial coordinates
    _x = p(0);
    _y = p(1);
    _z = p(2);
    
    // Evaluate time
    _t = t;

    // Evaluate Parameters (in case they changed since the constructor)
    for (auto it=_params.begin(); it!=_params.end(); it++) {
        _param_vars[it->first] = (*(it->second))();
    }

    uint n_components = _expression.size();
    // The symbol table is the same for all expressions so just grab it from any one of them
    exprtk::symbol_table<Real>& symbol_table = _dexpression.begin()->second[0].get_symbol_table();

    if (symbol_table.is_variable(f)) {
        // If the parameter we are taking the derivative with respect to (f) 
        // is a varaible, calculate the derivative.
        if (_dexpression.count(f)) {
            // Use the user-supplied analytical derivative if available
            for (uint i=0; i<n_components; i++) {
                v[i] = _dexpression.at(f)[i].value();
            }
        }
        else {
            // Use exprtk's finite difference approximation if analytical
            // derivative was not available.
            for (uint i=0; i<n_components; i++) {
                v[i] = exprtk::derivative(_expression[i], f, h);
            }
        }
    }
    else {
        // If the parameter we are taking the derivative with respect to (f)
        // is not a varaible, then the derivative is zero.
        for (uint i=0; i<n_components; i++) {
            v[i] = 0.0;
        }
    }
}


const std::vector<std::string> MAST::ExpressionVectorFieldFunction::get_expression() const {
    return _expr;
}