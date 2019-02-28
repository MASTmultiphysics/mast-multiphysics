/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2019  Manav Bhatia
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
#include "base/system_initialization.h"
#include "base/nonlinear_system.h"
#include "numerics/utility.h"
#include "base/field_function_base.h"

// libMesh includes
#include "libmesh/function_base.h"




MAST::SystemInitialization::SystemInitialization (MAST::NonlinearSystem& sys,
                                                  const std::string& prefix):
_system(sys),
_prefix(prefix) {

    // initialize the point locator for this mesh
    sys.system().get_mesh().sub_point_locator();
}



MAST::SystemInitialization::~SystemInitialization()
{ }




unsigned int
MAST::SystemInitialization::n_vars() const {
    return _system.n_vars();
}




const libMesh::FEType&
MAST::SystemInitialization::fetype(unsigned int i) const {
    
    return _system.variable_type(i);
}



void
MAST::SystemInitialization::initialize_solution(const RealVectorX& sol) {
    
    // make sure that the dimension of the sol vector matches the dimension
    // specified for this system
    libmesh_assert_equal_to(sol.size(), _vars.size());
    
    // now create a function and use it for initialization
    class SolutionFunction:
    public libMesh::FunctionBase<Real> {
    public:
        SolutionFunction(const RealVectorX& s):
        libMesh::FunctionBase<Real>() {
            _sol.resize((unsigned int)s.size());
            for (unsigned int i=0; i<s.size(); i++) _sol(i) = s(i);
        }
        
        
        SolutionFunction(const DenseRealVector& sol):
        libMesh::FunctionBase<Real>(),
        _sol(sol) { }
        
        virtual std::unique_ptr<libMesh::FunctionBase<Real> > clone () const {
            libMesh::FunctionBase<Real> *rval = new SolutionFunction(_sol);
            return std::unique_ptr<libMesh::FunctionBase<Real> >(rval);
        }
        
        // this should not get called
        virtual Real operator()
        (const libMesh::Point& p, const Real time) {libmesh_assert(false);}
        
        virtual void
        operator() (const libMesh::Point& p,
                    const Real time,
                    libMesh::DenseVector<Real>& output) {
            output = _sol;
        }
    protected:
        DenseRealVector _sol;
    };
    
    SolutionFunction sol_func(sol);
    
    _system.project_solution(&sol_func);
}



void
MAST::SystemInitialization::
initialize_solution(const MAST::FieldFunction<RealVectorX>& sol) {
    
    // now create a function and use it for initialization
    class SolutionFunction:
    public libMesh::FunctionBase<Real> {
    public:
        SolutionFunction(const unsigned int n_vars,
                         const MAST::FieldFunction<RealVectorX>& s):
        libMesh::FunctionBase<Real>(),
        _nvars(n_vars),
        _sol(s) { }
        
        virtual std::unique_ptr<libMesh::FunctionBase<Real> > clone () const {
            libMesh::FunctionBase<Real> *rval = new SolutionFunction(_nvars, _sol);
            return std::unique_ptr<libMesh::FunctionBase<Real> >(rval);
        }
        
        // this should not get called
        virtual Real operator()
        (const libMesh::Point& p, const Real time) {libmesh_assert(false);}
        
        virtual void
        operator() (const libMesh::Point& p,
                    const Real time,
                    libMesh::DenseVector<Real>& output) {
            
            RealVectorX v = RealVectorX::Zero(_nvars);
            _sol(p, time, v);
            MAST::copy(output, v);
        }
    protected:
        const unsigned int _nvars;
        const MAST::FieldFunction<RealVectorX>& _sol;
    };
    
    SolutionFunction sol_func(_vars.size(), sol);
    
    _system.project_solution(&sol_func);
}


