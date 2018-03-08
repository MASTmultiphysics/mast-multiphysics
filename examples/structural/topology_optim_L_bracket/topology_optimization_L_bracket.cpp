/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2018  Manav Bhatia
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
#include "examples/structural/topology_optim_L_bracket/topology_optimization_L_bracket.h"
#include "examples/base/input_wrapper.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "base/mesh_field_function.h"
#include "base/nonlinear_system.h"
#include "base/transient_assembly.h"
#include "base/boundary_condition_base.h"
#include "base/physics_discipline_base.h"
#include "boundary_condition/dirichlet_boundary_condition.h"

// libMesh includes
#include "libmesh/serial_mesh.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/elem.h"

MAST::Examples::TopologyOptimizationLevelSetLBracket::
TopologyOptimizationLevelSetLBracket(const libMesh::Parallel::Communicator& comm_in):
MAST::Examples::TopologyOptimizationLevelSet2D(comm_in) {
    
}



MAST::Examples::TopologyOptimizationLevelSetLBracket::
~TopologyOptimizationLevelSetLBracket() { }




void
MAST::Examples::TopologyOptimizationLevelSetLBracket::_init_mesh() {
    

    MAST::Examples::TopologyOptimizationLevelSet2D::_init_mesh();
    _delete_elems_from_mesh(*_mesh);
    _delete_elems_from_mesh(*_level_set_mesh);
}



void
MAST::Examples::TopologyOptimizationLevelSetLBracket::_init_dirichlet_conditions() {
 
    this->_init_boundary_dirichlet_constraint(0, "bottom_constraint");
    _discipline->init_system_dirichlet_bc(*_sys);
    _init_indicator_system_dirichlet_conditions();
}


namespace MAST {
    namespace Examples {
        class BracketLoad:
        public MAST::FieldFunction<Real> {
        public:
            BracketLoad(const std::string& nm, Real p, Real l1, Real fraction):
            MAST::FieldFunction<Real>(nm), _p(p), _l1(l1), _frac(fraction) { }
            virtual ~BracketLoad() {}
            virtual void operator() (const libMesh::Point& p, const Real t, Real& v) const {
                if (fabs(p(0) >= _l1*(1.-_frac))) v = _p;
                else v = 0.;
            }
            virtual void derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t, Real& v) const {
                v = 0.;
            }
        protected:
            Real _p, _l1, _frac;
        };
    }
}


void
MAST::Examples::TopologyOptimizationLevelSetLBracket::_init_loads() {

    Real
    length  = (*_input)(_prefix+"length", "length of domain along x-axis", 0.3),
    frac    = (*_input)(_prefix+"load_length_fraction", "fraction of boundary length on which pressure will act", 0.125),
    p_val   =  (*_input)(_prefix+"pressure", "pressure on side of domain",   5.e7);
    
    MAST::Examples::BracketLoad
    *press_f         = new MAST::Examples::BracketLoad( "pressure", p_val, length, frac),
    *flux_f          = new MAST::Examples::BracketLoad("heat_flux", -2.e6, length, frac);
    
    // initialize the load
    MAST::BoundaryConditionBase
    *p_load          = new MAST::BoundaryConditionBase(MAST::SURFACE_PRESSURE),
    *f_load          = new MAST::BoundaryConditionBase(MAST::HEAT_FLUX);
    
    p_load->add(*press_f);
    _discipline->add_side_load(2, *p_load);
    
    f_load->add(*flux_f);
    _indicator_discipline->add_side_load(2, *f_load);
    
    this->register_field_function(*press_f);
    this->register_field_function(*flux_f);
    this->register_loading(*p_load);
    this->register_loading(*f_load);
    
    _init_temperature_load();
}



void
MAST::Examples::TopologyOptimizationLevelSetLBracket::_init_phi_dvs() {
    
    libmesh_assert(_initialized);
    // this assumes that level set is defined using lagrange shape functions
    libmesh_assert_equal_to(_level_set_fetype.family, libMesh::LAGRANGE);
    
    Real
    tol     = 1.e-6,
    length  = (*_input)(_prefix+"length", "length of domain along x-axis", 0.3),
    height  = (*_input)(_prefix+"height", "length of domain along y-axis", 0.3),
    l_frac  = (*_input)(_prefix+"length_fraction", "fraction of length along x-axis that is in the bracket", 0.4),
    w_frac  = (*_input)(_prefix+ "width_fraction", "fraction of length along y-axis that is in the bracket", 0.4),
    x_lim   = length * l_frac,
    y_lim   =  height * (1.-w_frac),
    frac    = (*_input)(_prefix+"load_length_fraction", "fraction of boundary length on which pressure will act", 0.125);
    
    unsigned int
    dof_id  = 0;
    
    Real
    val     = 0.;
    
    // iterate over all the node values
    libMesh::MeshBase::const_node_iterator
    it  = _level_set_mesh->nodes_begin(),
    end = _level_set_mesh->nodes_end();
    
    // maximum number of dvs is the number of nodes on the level set function
    // mesh. We will evaluate the actual number of dvs
    _dv_params.reserve(_level_set_mesh->n_nodes());
    _n_vars = 0;

    for ( ; it!=end; it++) {
        
        const libMesh::Node& n = **it;

        dof_id                     = n.dof_number(0, 0, 0);

        if ((std::fabs(n(1)-height) > tol) ||
            (n(0) < length*(1.-frac))) {
            
            std::ostringstream oss;
            oss << "dv_" << _n_vars;
            val                        = _level_set_sys->solution->el(dof_id);
            
            // on the boundary, set everything to be zero, so that there
            // is always a boundary there that the optimizer can move
            if (n(0) < tol                     ||  // left boundary
                std::fabs(n(0) - length) < tol ||  // right boundary
                std::fabs(n(1) - height) < tol ||  // top boundary
                (n(0) >= x_lim && n(1) <= y_lim)) {
                
                _level_set_sys->solution->set(dof_id, 0.);
                val = 0.;
            }
            
            _dv_params.push_back(std::pair<unsigned int, MAST::Parameter*>());
            _dv_params[_n_vars].first  = dof_id;
            _dv_params[_n_vars].second = new MAST::Parameter(oss.str(), val);
            _dv_params[_n_vars].second->set_as_topology_parameter(true);
            
            _n_vars++;
        }
        else {
            // set value at the constrained points to a small positive number
            // material here
            _level_set_sys->solution->set(dof_id, 0.01);
        }
    }
    
    _level_set_sys->solution->close();
}


void
MAST::Examples::TopologyOptimizationLevelSetLBracket::_delete_elems_from_mesh(libMesh::MeshBase &mesh) {
    
    Real
    x       = -1.,
    y       = -1.,
    length  = (*_input)(_prefix+"length", "length of domain along x-axis", 0.3),
    width   = (*_input)(_prefix+ "width", "length of domain along y-axis", 0.3),
    l_frac  = (*_input)(_prefix+"length_fraction", "fraction of length along x-axis that is in the bracket", 0.4),
    w_frac  = (*_input)(_prefix+ "width_fraction", "fraction of length along y-axis that is in the bracket", 0.4),
    x_lim   = length * l_frac,
    y_lim   =  width * (1.-w_frac);
    
    // now, remove elements that are outside of the L-bracket domain
    libMesh::MeshBase::element_iterator
    e_it   = mesh.elements_begin(),
    e_end  = mesh.elements_end();
    
    for ( ; e_it!=e_end; e_it++) {
        
        libMesh::Elem* elem = *e_it;
        x = length;
        y = 0.;
        for (unsigned int i=0; i<elem->n_nodes(); i++) {
            const libMesh::Node& n = elem->node_ref(i);
            if (x > n(0)) x = n(0);
            if (y < n(1)) y = n(1);
        }
        
        // delete element if the lowest x,y locations are outside of the bracket
        // domain
        if (x >= x_lim && y<= y_lim)
            mesh.delete_elem(elem);
    }
    
    mesh.prepare_for_use();
}

