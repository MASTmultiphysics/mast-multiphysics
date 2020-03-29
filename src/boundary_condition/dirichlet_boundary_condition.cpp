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
#include "boundary_condition/dirichlet_boundary_condition.h"
#include "base/field_function_base.h"


// libMesh includes
#include "libmesh/zero_function.h"

void
MAST::DirichletBoundaryCondition::
init(const libMesh::boundary_id_type bid,
     const std::vector<unsigned int>& constrained_vars,
     MAST::FieldFunction<RealVectorX>* f_val,
     libMesh::VariableIndexing index,
     unsigned int n_sys_vars) {
    
    // should not have been initialized if this is called
    libmesh_assert(_dirichlet_boundary.get() == nullptr);
    
    std::set<libMesh::boundary_id_type> bid_set;
    bid_set.insert(bid);
    
    if (!f_val) {
        
        // if the function was not give, then assume it to be zero function
        std::unique_ptr<libMesh::FunctionBase<Real> >
        function(new libMesh::ZeroFunction<Real>);
                
        _dirichlet_boundary.reset(new libMesh::DirichletBoundary(bid_set,
                                                                 constrained_vars,
                                                                 function.get()));
    }
    else {
        
        // make sure that the number of variables is specified if indexing
        // is set to system_variable_order.
        if (index == libMesh::SYSTEM_VARIABLE_ORDER)
            libmesh_assert_greater(n_sys_vars, 0);
        
        // create a wrapper to the provided field function and pass it to the
        // DirichletBoundary object
        class FunctionWrapper: public libMesh::FunctionBase<Real> {
        public:
            FunctionWrapper(MAST::FieldFunction<RealVectorX>& f,
                            libMesh::VariableIndexing index,
                            unsigned int n_constrained_vars,
                            unsigned int n_sys_vars):
            libMesh::FunctionBase<Real>(),
            _f(f),
            _index(index),
            _n_constrained_vars(n_constrained_vars),
            _n_sys_vars(n_sys_vars) {
                
                if (_index == libMesh::SYSTEM_VARIABLE_ORDER)
                    libmesh_assert_greater(n_sys_vars, 0);
            }

            // copy constructor
            FunctionWrapper(const FunctionWrapper& f):
            libMesh::FunctionBase<Real>(),
            _f(f._f),
            _index(f._index),
            _n_constrained_vars(f._n_constrained_vars),
            _n_sys_vars(f._n_sys_vars) { }

            virtual ~FunctionWrapper() {}
            
            virtual std::unique_ptr<libMesh::FunctionBase<Real>> clone () const {
                
                std::unique_ptr<libMesh::FunctionBase<Real>> f;
                f.reset(new FunctionWrapper(*this));
                return f;
            }

            /**
             * \returns The scalar function value at coordinate \p p and time \p
             * time, which defaults to zero.
             *
             * Pure virtual, so you have to override it.
             */
            virtual Real operator() (const libMesh::Point & p,
                                     const Real time = 0.) {
                // should not get called
                libmesh_error();
            }

            virtual void operator() (const libMesh::Point & p,
                                     const Real time,
                                     libMesh::DenseVector<Real>& output) {
             
                RealVectorX v;
                _f(p, time, v);
                
                if (_index == libMesh::SYSTEM_VARIABLE_ORDER) {
                    libmesh_assert_equal_to(v.size(), _n_sys_vars);
                    output.resize(_n_sys_vars);
                    for (unsigned int i=0; i<_n_sys_vars; i++) output(i) = v(i);
                }
                else if (_index == libMesh::LOCAL_VARIABLE_ORDER) {
                 
                    libmesh_assert_equal_to(v.size(), _n_constrained_vars);
                    output.resize(_n_constrained_vars);
                    for (unsigned int i=0; i<_n_constrained_vars; i++) output(i) = v(i);
                }
            }

            
        protected:
            MAST::FieldFunction<RealVectorX> &_f;
            libMesh::VariableIndexing _index;
            unsigned int _n_constrained_vars;
            unsigned int _n_sys_vars;
        };
        
        std::unique_ptr<libMesh::FunctionBase<Real>>
        function(new FunctionWrapper(*f_val,
                                     index,
                                     constrained_vars.size(),
                                     n_sys_vars));
        
        _dirichlet_boundary.reset(new libMesh::DirichletBoundary(bid_set,
                                                                 constrained_vars,
                                                                 *function,
                                                                 index));
    }
}


