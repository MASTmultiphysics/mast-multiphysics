/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2016  Manav Bhatia
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
#include "fluid/conservative_fluid_system_initialization.h"


// libMesh includes
#include "libmesh/function_base.h"


MAST::ConservativeFluidSystemInitialization::
ConservativeFluidSystemInitialization(libMesh::System& sys,
                                      const std::string& prefix,
                                      const libMesh::FEType& fe_type,
                                      const unsigned int dim):
MAST::SystemInitialization(sys, prefix),
_dim(dim) {
    
    _vars.resize(dim+2);
    
    std::string nm = prefix + "_rho";
    _vars[0] = sys.add_variable(nm, fe_type);
    
    nm = prefix + "_rhoux";
    _vars[1] = sys.add_variable(nm, fe_type);
    
    if (dim > 1) {
        nm = prefix + "_rhouy";
        _vars[2] = sys.add_variable(nm, fe_type);
    }

    if (dim > 2) {
        nm = prefix + "_rhouz";
        _vars[3] = sys.add_variable(nm, fe_type);
    }

    
    nm = prefix + "_rhoe";
    _vars[dim+1] = sys.add_variable(nm, fe_type);
}



MAST::ConservativeFluidSystemInitialization::
~ConservativeFluidSystemInitialization() {
    
}




void
MAST::ConservativeFluidSystemInitialization::
initialize_solution(const RealVectorX& conservative_sol) {
    
    // make sure that the dimension of the sol vector matches the dimension
    // specified for this system
    libmesh_assert_equal_to(conservative_sol.size(), _vars.size());
    
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
        
        virtual libMesh::UniquePtr<FunctionBase<Real> > clone () const {
            FunctionBase<Real> *rval = new SolutionFunction(_sol);
            return libMesh::UniquePtr<FunctionBase<Real> >(rval);
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
    
    SolutionFunction sol_func(conservative_sol);

    _system.project_solution(&sol_func);
}



