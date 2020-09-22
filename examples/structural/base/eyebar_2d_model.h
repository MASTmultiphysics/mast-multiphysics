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

#ifndef __mast_topology_2d_eyebar_model__
#define __mast_topology_2d_eyebar_model__

// MAST includes
#include "base/mast_data_types.h"
#include "base/boundary_condition_base.h"
#include "base/field_function_base.h"
#include "base/physics_discipline_base.h"
#include "examples/base/input_wrapper.h"
#include "examples/structural/base/level_set_nucleation.h"
#include "examples/fluid/meshing/cylinder.h"
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
#include "libmesh/boundary_info.h"




namespace MAST {

// Forward declerations
class DisciplineBase;
class BoundaryConditionBase;
class FunctionBase;
class Parameter;
class Solid2DSectionElementPropertyCard;

namespace Examples {

struct Eyebar2DModel {
    
    using SectionPropertyCardType = MAST::Solid2DSectionElementPropertyCard;

    template <typename Opt>
    static Real reference_volume(Opt& opt);

    template <typename Opt>
    static void init_analysis_mesh(Opt& opt, libMesh::UnstructuredMesh& mesh);
    
    template <typename Opt>
    static void init_level_set_mesh(Opt& opt, libMesh::UnstructuredMesh& mesh);
    
    template <typename Opt>
    static void
    init_analysis_dirichlet_conditions(Opt& opt);
    
    template <typename Opt>
    static void
    init_indicator_dirichlet_conditions(Opt& opt);
    
    
    template <typename Opt>
    static MAST::BoundaryConditionBase&
    init_structural_shifted_boudnary_load(Opt& opt, unsigned int bid);

    template <typename Opt>
    static void
    init_structural_loads(Opt& opt);
    
    template <typename Opt>
    static void
    init_thermoelastic_loads(Opt& opt);

    template <typename Opt>
    static void
    init_indicator_loads(Opt& opt);
    
    template <typename Opt>
    static void
    init_level_set_dvs(Opt& opt);
    
    template <typename Opt>
    static void
    initialize_level_set_solution(Opt& opt);
    
    template <typename Opt>
    static void
    init_simp_dvs(Opt& opt);
    
    class EyebarLoad:
    public MAST::FieldFunction<Real> {
    public:
        EyebarLoad():
        MAST::FieldFunction<Real>("pressure") { }
        ~EyebarLoad() {}
        void operator() (const libMesh::Point& p, const Real t, Real& v) const {
            if (p(0) <= 0.) v = (-std::pow(p(1), 2) + std::pow(1.5, 2))*1.e6;
            else v = 0.;
        }
        void derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t, Real& v) const {
            v = 0.;
        }
    };
};

}
}


template <typename Opt>
Real
MAST::Examples::Eyebar2DModel::
reference_volume(Opt& opt) {
    
    return 16.*8. - acos(-1.) * 1.5*1.5;
}

    


template <typename Opt>
void
MAST::Examples::Eyebar2DModel::
init_analysis_mesh(Opt& opt, libMesh::UnstructuredMesh& mesh) {
    
    //
    // identify the element type from the input file or from the order
    // of the element
    //
    unsigned int
    n_radial_divs  = opt._input("n_radial_divs", "number of elements along radial direction", 20),
    n_quarter_divs = opt._input("n_quarter_divs", "number of elements along height", 20);
    
    Real
    //length   = 16.,
    height   = 8.,
    radius   = 1.5,
    h_ratio  = opt._input("h_ratio", "ratio of radial element size at cylinder and at edge", 2);
    
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
    
    MAST::Examples::CylinderMesh2D cylinder;
    cylinder.mesh(radius, height/2.,
                  n_radial_divs, n_quarter_divs, h_ratio,
                  mesh, e_type,
                  true, height, n_quarter_divs*2);
    
    //
    // add the boundary ids for Dirichlet conditions
    //
    libMesh::MeshBase::const_element_iterator
    e_it   = mesh.elements_begin(),
    e_end  = mesh.elements_end();
    
    Real
    tol  = radius * 1.e-8;
    
    for (; e_it != e_end; e_it++) {
        
        libMesh::Elem* elem = *e_it;
        
        std::unique_ptr<libMesh::Elem> edge(elem->side_ptr(1));
        libMesh::Point p = edge->centroid();
        
        if (std::fabs(p(0)-height*1.5) < tol &&
            std::fabs(p(1)) <= 1.) // on the right edge
            mesh.boundary_info->add_side(elem, 1, 0);
        
        // check for the circumference of the circle where load will be
        // applied
        edge.reset(elem->side_ptr(3).release());
        p = edge->centroid();
        
        if ((std::fabs(p.norm()-radius) < 1.e-2) &&
            p(0) < 0.) // left semi-circle
            mesh.boundary_info->add_side(elem, 3, 5);
    }
    
    mesh.boundary_info->sideset_name(0) = "dirichlet";
    mesh.boundary_info->sideset_name(5) = "load";
    
    
}


template <typename Opt>
void
MAST::Examples::Eyebar2DModel::
init_level_set_mesh(Opt& opt, libMesh::UnstructuredMesh& mesh) {
    
    Real
    //length   = 16.,
    height   = 8.,
    radius   = 1.5,
    h_ratio  = opt._input("h_ratio", "ratio of radial element size at cylinder and at edge", 2);
    
    unsigned int
    n_radial_divs  = opt._input("level_set_n_radial_divs", "number of elements along radial direction", 10),
    n_quarter_divs = opt._input("level_set_n_quarter_divs", "number of elements along height", 10);
    
    libMesh::ElemType
    e_type  = libMesh::QUAD4;
    
    MAST::Examples::CylinderMesh2D cylinder;
    cylinder.mesh(radius, height/2,
                  n_radial_divs, n_quarter_divs, h_ratio,
                  mesh, e_type,
                  true, height, n_quarter_divs*2);
}



template <typename Opt>
void
MAST::Examples::Eyebar2DModel::init_analysis_dirichlet_conditions(Opt& opt) {
    
    MAST::DirichletBoundaryCondition
    *dirichlet  = new MAST::DirichletBoundaryCondition;   // right boundary
    dirichlet->init(0, opt._sys_init->vars());
    opt._discipline->add_dirichlet_bc(0,  *dirichlet);
    opt._boundary_conditions.insert(dirichlet);
    
    opt._discipline->init_system_dirichlet_bc(*opt._sys);
    opt._dirichlet_bc_ids.insert(0);
}



template <typename Opt>
void
MAST::Examples::Eyebar2DModel::init_indicator_dirichlet_conditions(Opt& opt) {
    
    MAST::DirichletBoundaryCondition
    *dirichlet  = new MAST::DirichletBoundaryCondition;   // right boundary
    dirichlet->init(0, opt._indicator_sys_init->vars());
    opt._indicator_discipline->add_dirichlet_bc(0,  *dirichlet);
    opt._boundary_conditions.insert(dirichlet);
    
    opt._indicator_discipline->init_system_dirichlet_bc(*opt._indicator_sys);
}




template <typename Opt>
MAST::BoundaryConditionBase&
MAST::Examples::Eyebar2DModel::init_structural_shifted_boudnary_load(Opt& opt,
                                                                     unsigned int bid) {

    class ZeroTraction: public MAST::FieldFunction<RealVectorX> {
    public:
        ZeroTraction(): MAST::FieldFunction<RealVectorX>("traction") {}
        virtual ~ZeroTraction() {}
        virtual void operator() (const libMesh::Point& pt, const Real t, RealVectorX& v) const
        {v.setZero(3);}
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
MAST::Examples::Eyebar2DModel::init_structural_loads(Opt& opt) {
    
    EyebarLoad
    *press_f         = new EyebarLoad();
    
    // initialize the load
    MAST::BoundaryConditionBase
    *p_load          = new MAST::BoundaryConditionBase(MAST::SURFACE_PRESSURE);
    
    p_load->add(*press_f);
    opt._discipline->add_side_load(5, *p_load);
    opt._boundary_conditions.insert(p_load);

    opt._field_functions.insert(press_f);
}


template <typename Opt>
void
MAST::Examples::Eyebar2DModel::init_thermoelastic_loads(Opt& opt) {
    
    Real
    T_val      = opt._input(    "temperature",           "temperature for thermoelastic load",   0.),
    T_ref_val  = opt._input("ref_temperature", "reference temperature for thermoelastic load",   0.);

    MAST::Parameter
    *temperature = new MAST::Parameter("T",     T_val),
    *ref_temp    = new MAST::Parameter("T_ref", T_ref_val);

    MAST::ConstantFieldFunction
    *temperature_f = new MAST::ConstantFieldFunction(    "temperature", *temperature),
    *ref_temp_f    = new MAST::ConstantFieldFunction("ref_temperature",    *ref_temp);
    
    // initialize the load
    MAST::BoundaryConditionBase
    *t_load          = new MAST::BoundaryConditionBase(MAST::TEMPERATURE);

    t_load->add(*temperature_f);
    t_load->add(*ref_temp_f);
    opt._discipline->add_volume_load(0, *t_load);

    opt._boundary_conditions.insert(t_load);

    opt._parameters[temperature->name()] = temperature;
    opt._parameters[ref_temp->name()]    = ref_temp;
    opt._field_functions.insert(temperature_f);
    opt._field_functions.insert(ref_temp_f);
}




template <typename Opt>
void
MAST::Examples::Eyebar2DModel::init_indicator_loads(Opt& opt) {
    
}



template <typename Opt>
void
MAST::Examples::Eyebar2DModel::init_level_set_dvs(Opt& opt) {
    
    libmesh_assert(opt._initialized);
    //
    // this assumes that level set is defined using lagrange shape functions
    //
    libmesh_assert_equal_to(opt._level_set_fetype.family, libMesh::LAGRANGE);
    
    Real
    tol           = 1.e-6,
    height        = 8.,
    filterradius = opt._input("filterradius", "radius of geometric filter for level set field", 0.015);
    
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
        
        if (((n.norm() <= 1.5+filterradius) && n(0) <= 0.) ||  // circle
            (std::fabs(n(0)-height*1.5) < filterradius &&  // right edge
             std::fabs(n(1)) <= 1.+filterradius)) { // dirichlet constraint
            
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
            if (std::fabs(n(0)+height*0.5) < tol    ||  // left boundary
                std::fabs(n(1)-height*0.5) < tol    ||  // top boundary
                std::fabs(n(1)+height*0.5) < tol    ||  // bottom boundary
                std::fabs(n(0)-height*1.5) < tol) {     // right boundary
                
                if (dof_id >= opt._level_set_sys->solution->first_local_index() &&
                    dof_id <  opt._level_set_sys->solution->last_local_index())
                    opt._level_set_sys->solution->set(dof_id, -1.);
                val = -1.;
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
MAST::Examples::Eyebar2DModel::init_simp_dvs(Opt& opt) {
    
    libmesh_assert(opt._initialized);
    
    //
    // this assumes that density variable has a constant value per element
    //
    libmesh_assert_equal_to(opt._density_fetype.family, libMesh::LAGRANGE);
    
    Real
    tol           = 1.e-6,
    height        = 8.,
    filterradius  = opt._input("filterradius", "radius of geometric filter for level set field", 0.015);
    
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
    // iterate over all the element values
    libMesh::MeshBase::const_node_iterator
    it  = opt._mesh->nodes_begin(),
    end = opt._mesh->nodes_end();
    
    //
    // maximum number of dvs is the number of nodes on the level set function
    // mesh. We will evaluate the actual number of dvs
    //
    opt._dv_params.reserve(opt._mesh->n_elem());
    n_vars = 0;
    
    for ( ; it!=end; it++) {
        
        const libMesh::Node& n = **it;
        
        dof_id                     = n.dof_number(sys_num, 0, 0);
        
        
        
        if (((n.norm() <= 1.5+filterradius) && n(0) <= 0.) ||  // circle
            (std::fabs(n(0)-height*1.5) < filterradius &&  // right edge
             std::fabs(n(1)) <= 1.+filterradius)) { // dirichlet constraint
            
            //
            // set value at the constrained points to material
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
            if (std::fabs(n(0)+height*0.5) < tol    ||  // left boundary
                std::fabs(n(1)-height*0.5) < tol    ||  // top boundary
                std::fabs(n(1)+height*0.5) < tol    ||  // bottom boundary
                std::fabs(n(0)-height*1.5) < tol) {     // right boundary
                
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
MAST::Examples::Eyebar2DModel::initialize_level_set_solution(Opt& opt) {
    
    Real
    height = 8.,
    length = 16.;
    
    unsigned int
    nx_h    = opt._input("initial_level_set_n_holes_in_x",
                         "number of holes along x-direction for initial level-set field", 6),
    ny_h    = opt._input("initial_level_set_n_holes_in_y",
                         "number of holes along y-direction for initial level-set field", 6),
    nx_m    = opt._input("level_set_nx_divs", "number of elements of level-set mesh along x-axis", 10),
    ny_m    = opt._input("level_set_ny_divs", "number of elements of level-set mesh along y-axis", 10);
    
    MAST::Examples::LevelSetNucleationFunction2D
    phi(-0.5*height, -0.5*height, length, height, nx_m, ny_m, nx_h, ny_h);
    
    opt._level_set_sys_init->initialize_solution(phi);
}


#endif // __mast_topology_2d_eyebar_model__
