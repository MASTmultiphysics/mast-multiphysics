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
#include "base/expression_field_function.h"
#include "base/parameter.h"


MAST::ExpressionFieldFunction::
ExpressionFieldFunction(const std::string& nm,
                        const std::string expr,
                        const std::set<MAST::Parameter*> params):
MAST::FieldFunction<Real>(nm), 
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

    _expression.register_symbol_table(symbol_table);

    if (!(parser.compile(_expr, _expression)))
    {
        throw std::runtime_error("Error encoutered in exprtk when trying to compile expression.");
    }
}

MAST::ExpressionFieldFunction::
ExpressionFieldFunction(const std::string& nm,
                        const std::string expr,
                        const std::map<const std::string, const std::string> dexpr,
                        const std::set<MAST::Parameter*> params):
MAST::FieldFunction<Real>(nm), 
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

    _expression.register_symbol_table(symbol_table);
    if (!(parser.compile(_expr, _expression)))
    {
        throw std::runtime_error("Error encoutered in exprtk when trying to compile expression.");
    }

    std::map<const std::string, const std::string>::const_iterator it_expr;
    for (it_expr = _dexpr.begin(); it_expr != _dexpr.end(); it_expr++) {
        const std::string var = it_expr->first;
        const std::string dexpr_dvar = it_expr->second;
        exprtk::expression<Real> _expression;
        _dexpression.insert({var, _expression});
        _dexpression.at(var).register_symbol_table(symbol_table);
        if (!(parser.compile(dexpr_dvar, _dexpression.at(var))))
        {
            throw std::runtime_error("Error encoutered in exprtk when trying to compile dexpression.");
        }
    }
}


MAST::ExpressionFieldFunction::~ExpressionFieldFunction() {
    
}


void
MAST::ExpressionFieldFunction::operator() (const libMesh::Point& p, 
                                           const Real t, 
                                           Real& v) {
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

    v = _expression.value();
}


void
MAST::ExpressionFieldFunction::derivative (const MAST::FunctionBase& f,
                                           const libMesh::Point& p,
                                           const Real t,
                                           Real& v) {
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

    if (_expression.get_symbol_table().is_variable(f.name())) {
        if (_dexpression.count(f.name())) {
            v = _dexpression.at(f.name()).value();
        }
        else {
            v = exprtk::derivative(_expression, f.name()); // Finite difference approximation
        }
    }
    else {
        v = 0.0;
    }
}


void
MAST::ExpressionFieldFunction::derivative (const MAST::FunctionBase& f,
                                           const libMesh::Point& p,
                                           const Real t,
                                           Real& v,
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

    if (_expression.get_symbol_table().is_variable(f.name())) {
        if (_dexpression.count(f.name())) {
            v = _dexpression.at(f.name()).value();
        }
        else {
            v = exprtk::derivative(_expression, f.name(), h); // Finite difference approximation
        }
    }
    else {
        v = 0.0;
    }
}


void
MAST::ExpressionFieldFunction::derivative (const std::string& f,
                                           const libMesh::Point& p,
                                           const Real t,
                                           Real& v) {
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

    if (_expression.get_symbol_table().is_variable(f)) {
        if (_dexpression.count(f)) {
            v = _dexpression.at(f).value();
        }
        else {
            v = exprtk::derivative(_expression, f); // Finite difference approximation
        }
    }
    else {
        v = 0.0;
    }
}


void
MAST::ExpressionFieldFunction::derivative (const std::string& f,
                                           const libMesh::Point& p,
                                           const Real t,
                                           Real& v,
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

    if (_expression.get_symbol_table().is_variable(f)) {
        if (_dexpression.count(f)) {
            v = _dexpression.at(f).value();
        }
        else {
            v = exprtk::derivative(_expression, f, h); // Finite difference approximation
        }
    }
    else {
        v = 0.0;
    }
}


const std::string MAST::ExpressionFieldFunction::get_expression() const {
    return _expr;
}