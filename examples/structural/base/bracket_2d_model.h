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

#ifndef __mast_topology_2d_bracket_model__
#define __mast_topology_2d_bracket_model__

// MAST includes
#include "base/mast_data_types.h"
#include "examples/base/input_wrapper.h"
#include "examples/structural/base/level_set_nucleation.h"
#include "base/boundary_condition_base.h"
#include "base/field_function_base.h"
#include "base/physics_discipline_base.h"
#include "boundary_condition/dirichlet_boundary_condition.h"
#include "level_set/level_set_parameter.h"

// libMesh includes
#include "libmesh/system.h"
#include "libmesh/unstructured_mesh.h"
#include "libmesh/fe_type.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/elem.h"
#include "libmesh/node.h"



namespace MAST {

// Forward declerations
class DisciplineBase;
class BoundaryConditionBase;
class FunctionBase;
class Parameter;

namespace Examples {
struct Bracket2DModel {
    
    template <typename Opt>
    static Real reference_volume(Opt& opt);

    template <typename Opt>
    static void init_analysis_mesh(Opt& opt, libMesh::UnstructuredMesh& mesh);
    
    template <typename Opt>
    static void init_level_set_mesh(Opt& opt, libMesh::UnstructuredMesh& mesh);
    
    template <typename Opt>
    static void init_analysis_dirichlet_conditions(Opt& opt);
    
    template <typename Opt>
    static void init_indicator_dirichlet_conditions(Opt& opt);
    
    template <typename Opt>
    static void init_structural_loads(Opt& opt);

    template <typename Opt>
    static MAST::BoundaryConditionBase&
    init_structural_shifted_boudnary_load(Opt& opt, unsigned int bid);
    
    template <typename Opt>
    static void init_indicator_loads(Opt& opt);
    
    template <typename Opt>
    static void init_level_set_dvs(Opt& opt);
    
    template <typename Opt>
    static void initialize_level_set_solution(Opt& opt);
    
    template <typename Opt>
    static void init_simp_dvs(Opt& opt);
    
    template <typename Opt>
    static void _delete_elems_from_bracket_mesh(Opt& opt, libMesh::MeshBase &mesh);
    
    class BracketLoad:
    public MAST::FieldFunction<Real> {
    public:
        BracketLoad(const std::string& nm, Real p, Real l1, Real fraction):
        MAST::FieldFunction<Real>(nm), _p(p), _l1(l1), _frac(fraction) { }
        ~BracketLoad() {}
        void operator() (const libMesh::Point& p, const Real t, Real& v) const {
            if (fabs(p(0) >= _l1*(1.-_frac))) v = _p;
            else v = 0.;
        }
        void derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t, Real& v) const {
            v = 0.;
        }
    protected:
        Real _p, _l1, _frac;
    };
};

}
}


template <typename Opt>
Real
MAST::Examples::Bracket2DModel::
reference_volume(Opt& opt) {
    
    Real
    length  = opt._input("length", "length of domain along x-axis", 0.3),
    height  = opt._input("height", "length of domain along y-axis", 0.3);
    
    return length * height;
}

    
    
template <typename Opt>
void
MAST::Examples::Bracket2DModel::
init_analysis_mesh(Opt& opt,
                   libMesh::UnstructuredMesh& mesh) {
    
    Real
    length  = opt._input("length", "length of domain along x-axis", 0.3),
    height  = opt._input("height", "length of domain along y-axis", 0.3);
    
    unsigned int
    nx_divs = opt._input("nx_divs", "number of elements along x-axis", 20),
    ny_divs = opt._input("ny_divs", "number of elements along y-axis", 20);
    
    if (nx_divs%10 != 0 || ny_divs%10 != 0) libmesh_error();
    
    std::string
    t = opt._input("elem_type", "type of geometric element in the mesh", "quad4");
    
    libMesh::ElemType
    e_type = libMesh::Utility::string_to_enum<libMesh::ElemType>(t);
    
    //
    // if high order FE is used, libMesh requires atleast a second order
    // geometric element.
    //
    if (opt._fetype.order > 1 && e_type == libMesh::QUAD4)
        e_type = libMesh::QUAD9;
    else if (opt._fetype.order > 1 && e_type == libMesh::TRI3)
        e_type = libMesh::TRI6;
    
    //
    // initialize the mesh with one element
    //
    libMesh::MeshTools::Generation::build_square(mesh,
                                                 nx_divs, ny_divs,
                                                 0, length,
                                                 0, height,
                                                 e_type);
    
    _delete_elems_from_bracket_mesh(opt, mesh);
}


template <typename Opt>
void
MAST::Examples::Bracket2DModel::
init_level_set_mesh(Opt& opt,
                    libMesh::UnstructuredMesh& mesh) {
    
    Real
    length  = opt._input("length", "length of domain along x-axis", 0.3),
    height  = opt._input("height", "length of domain along y-axis", 0.3);
    
    unsigned int
    nx_divs = opt._input("level_set_nx_divs", "number of elements of level-set mesh along x-axis", 10),
    ny_divs = opt._input("level_set_ny_divs", "number of elements of level-set mesh along y-axis", 10);
    
    if (nx_divs%10 != 0 || ny_divs%10 != 0) libmesh_error();
    
    libMesh::ElemType
    e_type  = libMesh::QUAD4;
    
    // initialize the mesh with one element
    libMesh::MeshTools::Generation::build_square(mesh,
                                                 nx_divs, ny_divs,
                                                 0, length,
                                                 0, height,
                                                 e_type);
    
    _delete_elems_from_bracket_mesh(opt, mesh);
}



template <typename Opt>
void
MAST::Examples::Bracket2DModel::
init_analysis_dirichlet_conditions(Opt& opt) {
    
    MAST::DirichletBoundaryCondition
    *dirichlet  = new MAST::DirichletBoundaryCondition;   // bottom boundary
    dirichlet->init(0, opt._sys_init->vars());
    opt._discipline->add_dirichlet_bc(0,  *dirichlet);
    opt._boundary_conditions.insert(dirichlet);
    
    opt._discipline->init_system_dirichlet_bc(*opt._sys);
}



template <typename Opt>
void
MAST::Examples::Bracket2DModel::
init_indicator_dirichlet_conditions(Opt& opt) {
    
    MAST::DirichletBoundaryCondition
    *dirichlet  = new MAST::DirichletBoundaryCondition;   // bottom boundary
    dirichlet->init(0, opt._indicator_sys_init->vars());
    opt._indicator_discipline->add_dirichlet_bc(0,  *dirichlet);
    opt._boundary_conditions.insert(dirichlet);
    
    opt._indicator_discipline->init_system_dirichlet_bc(*opt._indicator_sys);
    opt._dirichlet_bc_ids.insert(0);
}



template <typename Opt>
void
MAST::Examples::Bracket2DModel::init_structural_loads(Opt& opt) {
    
    Real
    length  = opt._input("length", "length of domain along x-axis", 0.3),
    frac    = opt._input("loadlength_fraction", "fraction of boundary length on which pressure will act", 0.125),
    p_val   = opt._input("pressure", "pressure on side of domain",   5.e7);
    
    BracketLoad
    *press_f         = new BracketLoad( "pressure", p_val, length, frac);
    
    //
    // initialize the load
    //
    MAST::BoundaryConditionBase
    *p_load          = new MAST::BoundaryConditionBase(MAST::SURFACE_PRESSURE);
    
    p_load->add(*press_f);
    opt._discipline->add_side_load(5, *p_load);
    opt._boundary_conditions.insert(p_load);

    opt._field_functions.insert(press_f);
}


template <typename Opt>
MAST::BoundaryConditionBase&
MAST::Examples::Bracket2DModel::init_structural_shifted_boudnary_load(Opt& opt,
                                                                      unsigned int bid) {

    class ZeroTraction: public MAST::FieldFunction<RealVectorX> {
    public:
        ZeroTraction(): MAST::FieldFunction<RealVectorX>("traction") {}
        virtual ~ZeroTraction() {}
        virtual void operator() (const libMesh::Point& pt, const Real t, RealVectorX& v) const        {v.setZero(3);}
        virtual void derivative(const MAST::FunctionBase& f, const libMesh::Point& pt, const Real t, RealVectorX& v) const
        {v.setZero(3);}
    };
    
    ZeroTraction
    *trac_f = new ZeroTraction;
    
    MAST::BoundaryConditionBase
    *load          = new MAST::BoundaryConditionBase(MAST::SURFACE_TRACTION_SHIFTED_BOUNDARY);
    
    load->add(*opt._level_set_vel);
    load->add(*trac_f);
    opt._discipline->add_side_load(bid, *load);
    opt._boundary_conditions.insert(load);

    opt._field_functions.insert(trac_f);
    return *load;
}


template <typename Opt>
void
MAST::Examples::Bracket2DModel::init_indicator_loads(Opt& opt) {
    
    Real
    length  = opt._input("length", "length of domain along x-axis", 0.3),
    frac    = opt._input("loadlength_fraction", "fraction of boundary length on which pressure will act", 0.125);
    
    BracketLoad
    *flux_f          = new BracketLoad("heat_flux", -2.e6, length, frac);
    
    //
    // initialize the load
    //
    MAST::BoundaryConditionBase
    *f_load          = new MAST::BoundaryConditionBase(MAST::HEAT_FLUX);
    
    f_load->add(*flux_f);
    opt._indicator_discipline->add_side_load(5, *f_load);
    opt._boundary_conditions.insert(f_load);

    opt._field_functions.insert(flux_f);
}



template <typename Opt>
void
MAST::Examples::Bracket2DModel::init_level_set_dvs(Opt& opt) {
    
    libmesh_assert(opt._initialized);
    //
    // this assumes that level set is defined using lagrange shape functions
    //
    libmesh_assert_equal_to(opt._level_set_fetype.family, libMesh::LAGRANGE);
    
    Real
    tol           = 1.e-12,
    l_frac        = 0.4,//_input("length_fraction", "fraction of length along x-axis that is in the bracket", 0.4),
    h_frac        = 0.4,//_input( "height_fraction", "fraction of length along y-axis that is in the bracket", 0.4),
    length        = opt._input("length", "length of domain along x-axis", 0.3),
    height        = opt._input("height", "length of domain along y-axis", 0.3),
    x_lim         = length * l_frac,
    y_lim         = height * (1.-h_frac),
    frac          = opt._input("loadlength_fraction", "fraction of boundary length on which pressure will act", 0.125),
    filter_radius = opt._input("filter_radius", "radius of geometric filter for level set field", 0.015);
    
    unsigned int
    dof_id  = 0,
    n_vars  = 0;
    
    Real
    val     = 0.;
    
    //
    // all ranks will have DVs defined for all variables. So, we should be
    // operating on a replicated mesh
    //
    libmesh_assert(opt._level_set_mesh->is_replicated());
    
    std::vector<Real> local_phi(opt._level_set_sys->solution->size());
    opt._level_set_sys->solution->localize(local_phi);
    
    //
    // iterate over all the node values
    //
    libMesh::MeshBase::const_node_iterator
    it  = opt._level_set_mesh->nodes_begin(),
    end = opt._level_set_mesh->nodes_end();
    
    //
    // maximum number of dvs is the number of nodes on the level set function
    // mesh. We will evaluate the actual number of dvs
    //
    opt._dv_params.reserve(opt._level_set_mesh->n_nodes());
    n_vars = 0;
    
    for ( ; it!=end; it++) {
        
        const libMesh::Node& n = **it;
        
        dof_id                     = n.dof_number(0, 0, 0);
        
        if ((n(1)-filter_radius) <= y_lim &&
            (n(0)+filter_radius) >= length*(1.-frac)) {
            
            //
            // set value at the constrained points to a small positive number
            // material here
            //
            if (dof_id >= opt._level_set_sys->solution->first_local_index() &&
                dof_id <  opt._level_set_sys->solution->last_local_index())
                opt._level_set_sys->solution->set(dof_id, 1.e0);
        }
        else {
            
            std::ostringstream oss;
            oss << "dv_" << n_vars;
            val = local_phi[dof_id];
            
            //
            // on the boundary, set everything to be zero, so that there
            // is always a boundary there that the optimizer can move
            //
            if (n(0) < tol                     ||  // left boundary
                std::fabs(n(0) - length) < tol ||  // right boundary
                std::fabs(n(1) - height) < tol ||  // top boundary
                (n(0) >= x_lim && n(1) <= y_lim)) {
                
                if (dof_id >= opt._level_set_sys->solution->first_local_index() &&
                    dof_id <  opt._level_set_sys->solution->last_local_index())
                    opt._level_set_sys->solution->set(dof_id, -1.0);
                val = -1.0;
            }
            
            opt._dv_params.push_back(std::pair<unsigned int, MAST::Parameter*>());
            opt._dv_params[n_vars].first  = dof_id;
            opt._dv_params[n_vars].second = new MAST::LevelSetParameter(oss.str(), val, &n);
            opt._dv_params[n_vars].second->set_as_topology_parameter(true);
            opt._dv_dof_ids.insert(dof_id);
            
            n_vars++;
        }
    }
    
    opt.set_n_vars(n_vars);
    
    opt._level_set_sys->solution->close();
}


template <typename Opt>
void
MAST::Examples::Bracket2DModel::init_simp_dvs(Opt& opt) {
    
    libmesh_assert(opt._initialized);
    
    //
    // this assumes that density variable has a constant value per element
    //
    libmesh_assert_equal_to(opt._density_fetype.family, libMesh::LAGRANGE);
    
    Real
    tol           = 1.e-12,
    l_frac        = 0.4,//_input("length_fraction", "fraction of length along x-axis that is in the bracket", 0.4),
    h_frac        = 0.4,//_input( "height_fraction", "fraction of length along y-axis that is in the bracket", 0.4),
    length        = opt._input("length", "length of domain along x-axis", 0.3),
    height        = opt._input("height", "length of domain along y-axis", 0.3),
    x_lim         = length * l_frac,
    y_lim         = height * (1.-h_frac),
    frac          = opt._input("loadlength_fraction", "fraction of boundary length on which pressure will act", 0.125),
    filter_radius = opt._input("filter_radius", "radius of geometric filter for level set field", 0.015);
    
    unsigned int
    sys_num = opt._density_sys->number(),
    dof_id  = 0,
    n_vars  = 0;
    
    Real
    val     = 0.;
    
    //
    // all ranks will have DVs defined for all variables. So, we should be
    // operating on a replicated mesh
    //
    libmesh_assert(opt._mesh->is_replicated());
    
    std::vector<Real> local_phi(opt._density_sys->solution->size());
    opt._density_sys->solution->localize(local_phi);
    
    // iterate over all the element values
    libMesh::MeshBase::const_node_iterator
    it  = opt._mesh->nodes_begin(),
    end = opt._mesh->nodes_end();
    
    //
    // maximum number of dvs is the number of nodes on the level set function
    // mesh. We will evaluate the actual number of dvs
    //
    opt._dv_params.reserve(opt._mesh->n_elem());
    
    for ( ; it!=end; it++) {
        
        const libMesh::Node& n = **it;
        
        dof_id                     = n.dof_number(sys_num, 0, 0);
        
        if ((n(1)-filter_radius) <= y_lim && (n(0)+filter_radius) >= length*(1.-frac)) {
            
            //
            // set value at the constrained points to a small positive number
            // material here
            //
            if (dof_id >= opt._density_sys->solution->first_local_index() &&
                dof_id <  opt._density_sys->solution->last_local_index())
                opt._density_sys->solution->set(dof_id, 1.e0);
        }
        else {
            
            std::ostringstream oss;
            oss << "dv_" << n_vars;
            val = local_phi[dof_id];
            
            //
            // on the boundary, set everything to be zero, so that there
            // is always a boundary there that the optimizer can move
            //
            if (n(0) < tol                     ||  // left boundary
                std::fabs(n(0) - length) < tol ||  // right boundary
                std::fabs(n(1) - height) < tol ||  // top boundary
                (n(0) >= x_lim && n(1) <= y_lim)) {
                
                if (dof_id >= opt._density_sys->solution->first_local_index() &&
                    dof_id <  opt._density_sys->solution->last_local_index())
                    opt._density_sys->solution->set(dof_id, opt._rho_min);
                val = opt._rho_min;
            }
            
            opt._dv_params.push_back(std::pair<unsigned int, MAST::Parameter*>());
            opt._dv_params[n_vars].first  = dof_id;
            opt._dv_params[n_vars].second = new MAST::LevelSetParameter(oss.str(), val, &n);
            opt._dv_params[n_vars].second->set_as_topology_parameter(true);
            opt._dv_dof_ids.insert(dof_id);
            
            n_vars++;
        }
    }
    
    opt.set_n_vars(n_vars);
    opt._density_sys->solution->close();
    
}



template <typename Opt>
void
MAST::Examples::Bracket2DModel::
_delete_elems_from_bracket_mesh(Opt& opt,
                                libMesh::MeshBase &mesh) {
    
    Real
    tol     = 1.e-12,
    x       = -1.,
    y       = -1.,
    l_frac  = 0.4,
    w_frac  = 0.4,
    length  = opt._input("length", "length of domain along x-axis", 0.3),
    height  = opt._input("height", "length of domain along y-axis", 0.3),
    x_lim   = length * l_frac,
    y_lim   = height * (1.-w_frac);
    
    //
    // now, remove elements that are outside of the L-bracket domain
    //
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
        
        //
        // delete element if the lowest x,y locations are outside of the bracket
        // domain
        //
        if (x >= x_lim && y<= y_lim)
            mesh.delete_elem(elem);
    }
    
    mesh.prepare_for_use();
    
    //
    // add the two additional boundaries to the boundary info so that
    // we can apply loads on them
    //
    bool
    facing_right = false,
    facing_down  = false;
    
    e_it   = mesh.elements_begin();
    e_end  = mesh.elements_end();
    
    for ( ; e_it != e_end; e_it++) {
        
        libMesh::Elem* elem = *e_it;
        
        if (!elem->on_boundary()) continue;
        
        for (unsigned int i=0; i<elem->n_sides(); i++) {
            
            if (elem->neighbor_ptr(i)) continue;
            
            std::unique_ptr<libMesh::Elem> s(elem->side_ptr(i).release());
            
            const libMesh::Point p = s->centroid();
            
            facing_right = true;
            facing_down  = true;
            for (unsigned int j=0; j<s->n_nodes(); j++) {
                const libMesh::Node& n = s->node_ref(j);
                
                if (n(0) < x_lim ||  n(1) > y_lim) {
                    facing_right = false;
                    facing_down  = false;
                }
                else if (std::fabs(n(0) - p(0)) > tol)
                    facing_right = false;
                else if (std::fabs(n(1) - p(1)) > tol)
                    facing_down = false;
            }
            
            if (facing_right) mesh.boundary_info->add_side(elem, i, 4);
            if (facing_down) mesh.boundary_info->add_side(elem, i, 5);
        }
    }
    
    mesh.boundary_info->sideset_name(4) = "facing_right";
    mesh.boundary_info->sideset_name(5) = "facing_down";
}



template <typename Opt>
void
MAST::Examples::Bracket2DModel::initialize_level_set_solution(Opt& opt) {
    
    Real
    length  = opt._input("length", "length of domain along x-axis", 0.3),
    height  = opt._input("height", "length of domain along y-axis", 0.3);
    
    unsigned int
    nx_h    = opt._input("initial_level_set_n_holes_in_x",
                         "number of holes along x-direction for initial level-set field", 6),
    ny_h    = opt._input("initial_level_set_n_holes_in_y",
                         "number of holes along y-direction for initial level-set field", 6),
    nx_m    = opt._input("level_set_nx_divs", "number of elements of level-set mesh along x-axis", 10),
    ny_m    = opt._input("level_set_ny_divs", "number of elements of level-set mesh along y-axis", 10);
    
    MAST::Examples::LevelSetNucleationFunction
    phi(0., 0., length, height, nx_m, ny_m, nx_h, ny_h);
    
    opt._level_set_sys_init->initialize_solution(phi);
}


#endif // __mast_topology_2d_bracket_model__
