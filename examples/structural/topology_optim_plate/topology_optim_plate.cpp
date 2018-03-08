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
#include "examples/structural/topology_optim_plate/topology_optim_plate.h"
#include "examples/base/input_wrapper.h"
#include "boundary_condition/dirichlet_boundary_condition.h"
#include "base/physics_discipline_base.h"
#include "base/nonlinear_system.h"
#include "base/parameter.h"
#include "heat_conduction/heat_conduction_system_initialization.h"

// libMesh includes
#include "libmesh/node.h"
#include "libmesh/numeric_vector.h"


MAST::Examples::TopologyOptimizationLevelSetPlate::
TopologyOptimizationLevelSetPlate(const libMesh::Parallel::Communicator& comm_in):
MAST::Examples::TopologyOptimizationLevelSet2D(comm_in) {
    
}



MAST::Examples::TopologyOptimizationLevelSetPlate::~TopologyOptimizationLevelSetPlate() { }



void
MAST::Examples::TopologyOptimizationLevelSetPlate::_init_dirichlet_conditions() {
    
    MAST::Examples::StructuralExample2D::_init_dirichlet_conditions();
    
    ///////////////////////////////////////////////////////////////////////
    // initialize Dirichlet conditions for indicator system
    ///////////////////////////////////////////////////////////////////////
    MAST::DirichletBoundaryCondition
    *dirichlet  = new MAST::DirichletBoundaryCondition;   // right boundary
    dirichlet->init(1, _indicator_sys_init->vars());
    _indicator_discipline->add_dirichlet_bc(1,  *dirichlet);
    this->register_loading(*dirichlet);                   // register, so that it will be deleted later
    dirichlet   = new MAST::DirichletBoundaryCondition;   // left boundary
    dirichlet->init(3, _indicator_sys_init->vars());
    _indicator_discipline->add_dirichlet_bc(3,  *dirichlet);
    this->register_loading(*dirichlet);                   // register, so that it will be deleted later
    _indicator_discipline->init_system_dirichlet_bc(*_indicator_sys);
}



namespace MAST {
    namespace Examples {
        class PlateFluxLoadCenter:
        public MAST::FieldFunction<Real> {
        public:
            PlateFluxLoadCenter(const std::string& nm, Real p, Real l1, Real l2, Real fraction):
            MAST::FieldFunction<Real>(nm), _p(p), _l1(l1), _l2(l2), _frac(fraction) { }
            virtual ~PlateFluxLoadCenter() {}
            virtual void operator() (const libMesh::Point& p, const Real t, Real& v) const {
                Real d = pow(pow(p(0)-_l1*0.5, 2) + pow(p(1)-_l2*0.5, 2), 0.5);
                if (std::fabs(d - 0.25*_frac*(_l1+_l2)) <= 0.) v = _p;
                else v = 0.;
            }
            virtual void derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t, Real& v) const {
                v = 0.;
            }
        protected:
            Real _p, _l1, _l2, _frac;
        };
    }
}

void
MAST::Examples::TopologyOptimizationLevelSetPlate::_init_loads() {
    
    Real
    length  = (*_input)(_prefix+"length", "length of domain along x-axis", 0.3),
    height  = (*_input)(_prefix+"height", "length of domain along y-axis", 0.3),
    frac    = (*_input)(_prefix+"load_length_fraction", "fraction of length that identifies diameter over which pressure acts", 0.2),
    p_val   = (*_input)(_prefix+"pressure", "pressure on domain",   2.e4);
    
    MAST::Examples::PlateFluxLoadCenter
    *press_f         = new MAST::Examples::PlateFluxLoadCenter( "pressure", p_val, length, height, frac),
    *flux_f          = new MAST::Examples::PlateFluxLoadCenter("heat_flux", -2.e6, length, height, frac);
    
    // initialize the load
    MAST::BoundaryConditionBase
    *p_load          = new MAST::BoundaryConditionBase(MAST::SURFACE_PRESSURE),
    *f_load          = new MAST::BoundaryConditionBase(MAST::HEAT_FLUX);
    
    p_load->add(*press_f);
    _discipline->add_volume_load(0, *p_load);
    
    f_load->add(*flux_f);
    _indicator_discipline->add_volume_load(0, *f_load);
    
    this->register_field_function(*press_f);
    this->register_field_function(*flux_f);
    this->register_loading(*p_load);
    this->register_loading(*f_load);
    
    _init_temperature_load();
}



void
MAST::Examples::TopologyOptimizationLevelSetPlate::_init_phi_dvs() {
    
    libmesh_assert(_initialized);
    // this assumes that level set is defined using lagrange shape functions
    libmesh_assert_equal_to(_level_set_fetype.family, libMesh::LAGRANGE);
    
    Real
    tol     = 1.e-6,
    length  = (*_input)(_prefix+"length", "length of domain along x-axis", 0.3),
    height  = (*_input)(_prefix+"height", "length of domain along y-axis", 0.3),
    frac    = (*_input)(_prefix+"load_length_fraction", "fraction of boundary length on which pressure will act", 0.2);
    
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
        
        Real d = pow(pow(n(0)-length*0.5, 2) + pow(n(1)-length*0.5, 2), 0.5);
        dof_id                     = n.dof_number(0, 0, 0);

        // only if node is not on the upper edge
        if (/*n(0) > 0.     &&
            n(1) > 0.     &&
            n(0) < length &&
            n(1) < height &&*/
            std::fabs(d - 0.25*frac*(length+height)) >= 0.) {
         
            std::ostringstream oss;
            oss << "dv_" << _n_vars;
            val                        = _level_set_sys->solution->el(dof_id);
            
            // on the boundary, set everything to be zero, so that there
            // is always a boundary there
            /*if (n(0) < tol                     ||
                n(1) < tol                     ||
                std::fabs(n(0) - length) < tol ||
                std::fabs(n(1) - height) < tol) {
                
                _level_set_sys->solution->set(dof_id, 0.);
                val = 0.;
            }*/

            _dv_params.push_back(std::pair<unsigned int, MAST::Parameter*>());
            _dv_params[_n_vars].first  = dof_id;
            _dv_params[_n_vars].second = new MAST::Parameter(oss.str(), val);
            _dv_params[_n_vars].second->set_as_topology_parameter(true);
            
            _n_vars++;
        }
        else if (std::fabs(d - 0.25*frac*(length+height)) <= 0.) {
            // set value inside the circle to be 1, so that we always create a
            // material here
            _level_set_sys->solution->set(dof_id, 0.01);
        }
    }
    
    _level_set_sys->solution->close();
}



