/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2017  Manav Bhatia
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


// C++ includes
#include <iostream>


//MAST includes
#include "examples/structural/nastran_model_analysis/mast_interface.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "property_cards/isotropic_material_property_card.h"
#include "property_cards/solid_2d_section_element_property_card.h"
#include "base/nonlinear_system.h"
#include "elasticity/structural_discipline.h"
#include "elasticity/structural_system_initialization.h"


MAST::Model::Model(libMesh::LibMeshInit& init):
_mesh(new libMesh::SerialMesh(init.comm())) {

    MAST::Parameter *z = new MAST::Parameter("zero",  0.);
    _param_map["zero"] = z;
}



MAST::Model::~Model() {
    
    delete _mesh;

    {
        std::map<std::string, MAST::Parameter*>::iterator
        it  = _param_map.begin(),
        end = _param_map.end();
        
        for ( ; it != end; it++)
            delete it->second;
    }
    
    {
        std::map<std::string, MAST::FunctionBase*>::iterator
        it  = _function_map.begin(),
        end = _function_map.end();
        
        for ( ; it != end; it++)
            delete it->second;
    }
    
    
    {
        std::map<unsigned int, MAST::MaterialPropertyCardBase*>::iterator
        it  = _material_map.begin(),
        end = _material_map.end();
        
        for ( ; it != end; it++)
            delete it->second;
    }
    
    {
        std::map<unsigned int, MAST::ElementPropertyCardBase*>::iterator
        it  = _elem_property_map.begin(),
        end = _elem_property_map.end();
        
        for ( ; it != end; it++)
            delete it->second;
    }
    

    delete _discipline;
    delete _sys_init;
    delete _eq_sys;
}




void
MAST::Model::add_node(int i, double x, double y, double z) {

    // make sure this does not already exist
    libmesh_assert(!_node_id_map.count(i));
    
    // create the node
    libMesh::Node
    *n = _mesh->add_point(libMesh::Point(x, y, z));
    
    // add node to the map
    _node_id_map[i] = n;
}


void
MAST::Model::add_edge2(int e_id, int p_id, int n1, int n2) {

    // make sure this does not already exist
    libmesh_assert(!_elem_id_map.count(e_id));
    
    // create the element
    libMesh::Elem* e = libMesh::Elem::build(libMesh::EDGE2).release();
    e->subdomain_id() = p_id;
    
    // make sure that the nodes are contained in the map
    libmesh_assert(_node_id_map.count(n1));
    libmesh_assert(_node_id_map.count(n2));
    
    e->set_node(0) = _node_id_map[n1];
    e->set_node(1) = _node_id_map[n2];
    
    _mesh->add_elem(e);
    
    // add this elem to the map
    _elem_id_map[e_id] = e;
}



void
MAST::Model::add_quad4(int e_id, int p_id,
                       int n1, int n2,
                       int n3, int n4) {

    // make sure this does not already exist
    libmesh_assert(!_elem_id_map.count(e_id));
    
    // create the element
    libMesh::Elem* e = libMesh::Elem::build(libMesh::QUAD4).release();
    e->subdomain_id() = p_id;
    
    // make sure that the nodes are contained in the map
    libmesh_assert(_node_id_map.count(n1));
    libmesh_assert(_node_id_map.count(n2));
    libmesh_assert(_node_id_map.count(n3));
    libmesh_assert(_node_id_map.count(n4));
    
    e->set_node(0) = _node_id_map[n1];
    e->set_node(1) = _node_id_map[n2];
    e->set_node(2) = _node_id_map[n3];
    e->set_node(3) = _node_id_map[n4];
    
    _mesh->add_elem(e);
    
    // add this elem to the map
    _elem_id_map[e_id] = e;
}


void
MAST::Model::add_isotropic_material(int m_id, double E, double nu, double rho, double alpha) {
    
    // make sure this does not already exist
    libmesh_assert(!_material_map.count(m_id));
    
    // create the parameters for each value
    MAST::Parameter
    *E_p     = nullptr,
    *nu_p    = nullptr,
    *rho_p   = nullptr,
    *alpha_p = nullptr,
    *kappa_p = nullptr;

    MAST::ConstantFieldFunction
    *E_f     = nullptr,
    *nu_f    = nullptr,
    *rho_f   = nullptr,
    *alpha_f = nullptr,
    *kappa_f = nullptr;
    
    std::ostringstream oss;

    // create the parameters and functions for each variable
    oss << "Mat_" << m_id << "_E";
    E_p = new MAST::Parameter(oss.str(), E);
    E_f = new MAST::ConstantFieldFunction("E", *E_p);
    libmesh_assert(!_param_map.count(oss.str()));
    libmesh_assert(!_function_map.count(oss.str()));
    _param_map[oss.str()]    = E_p;
    _function_map[oss.str()] = E_f;
    
    oss.clear();
    oss << "Mat_" << m_id << "_nu";
    nu_p = new MAST::Parameter(oss.str(), nu);
    nu_f = new MAST::ConstantFieldFunction("nu", *nu_p);
    libmesh_assert(!_param_map.count(oss.str()));
    libmesh_assert(!_function_map.count(oss.str()));
    _param_map[oss.str()]    = nu_p;
    _function_map[oss.str()] = nu_f;

    oss.clear();
    oss << "Mat_" << m_id << "_rho";
    rho_p = new MAST::Parameter(oss.str(), rho);
    rho_f = new MAST::ConstantFieldFunction("rho", *rho_p);
    libmesh_assert(!_param_map.count(oss.str()));
    libmesh_assert(!_function_map.count(oss.str()));
    _param_map[oss.str()]    = rho_p;
    _function_map[oss.str()] = rho_f;

    oss.clear();
    oss << "Mat_" << m_id << "_alpha";
    alpha_p = new MAST::Parameter(oss.str(), alpha);
    alpha_f = new MAST::ConstantFieldFunction("alpha", *alpha_p);
    libmesh_assert(!_param_map.count(oss.str()));
    libmesh_assert(!_function_map.count(oss.str()));
    _param_map[oss.str()]    = alpha_p;
    _function_map[oss.str()] = alpha_f;

    oss.clear();
    oss << "Mat_" << m_id << "_kappa";
    kappa_p = new MAST::Parameter(oss.str(), 5./6.);
    kappa_f = new MAST::ConstantFieldFunction("kappa", *kappa_p);
    libmesh_assert(!_param_map.count(oss.str()));
    libmesh_assert(!_function_map.count(oss.str()));
    _param_map[oss.str()]    = kappa_p;
    _function_map[oss.str()] = kappa_f;
    
    // create the material property card
    MAST::IsotropicMaterialPropertyCard
    *mat = new MAST::IsotropicMaterialPropertyCard;
    mat->add(*E_f);
    mat->add(*nu_f);
    mat->add(*rho_f);
    mat->add(*alpha_f);
    mat->add(*kappa_f);
    _material_map[m_id] = mat;
    
}



void
MAST::Model::add_2d_section_property(int p_id, int m_id, double t) {
    
    // make sure this does not already exist
    libmesh_assert(!_elem_property_map.count(p_id));
    libmesh_assert(_material_map.count(m_id));
    
    // create the parameters for each value
    MAST::Parameter
    *h_p     = nullptr,
    *zero    = _param_map["zero"];
    
    MAST::ConstantFieldFunction
    *h_f     = nullptr,
    *off_f   = nullptr;
    
    std::ostringstream oss;
    
    // create the parameters and functions for each variable
    oss << "Prop_" << p_id << "_h";
    h_p = new MAST::Parameter(oss.str(), t);
    h_f = new MAST::ConstantFieldFunction("h", *h_p);
    libmesh_assert(!_param_map.count(oss.str()));
    libmesh_assert(!_function_map.count(oss.str()));
    _param_map[oss.str()]    = h_p;
    _function_map[oss.str()] = h_f;
    
    oss.clear();
    oss << "Prop_" << p_id << "_off";
    off_f = new MAST::ConstantFieldFunction("off", *zero);
    libmesh_assert(!_param_map.count(oss.str()));
    libmesh_assert(!_function_map.count(oss.str()));
    _function_map[oss.str()] = off_f;
    
    // create the material property card
    MAST::Solid2DSectionElementPropertyCard
    *prop = new MAST::Solid2DSectionElementPropertyCard;
    MAST::MaterialPropertyCardBase
    *mat  = _material_map[m_id];
    prop->add(*h_f);
    prop->add(*off_f);
    prop->set_material(*mat);
    
    _elem_property_map[p_id] = prop;

}


void
MAST::Model::initialize_after_mesh() {
    
    // prepare the mesh for use
    _mesh->prepare_for_use();
    
    _eq_sys      = new  libMesh::EquationSystems(*_mesh);
    
    // create the libmesh system
    _sys         = &(_eq_sys->add_system<MAST::NonlinearSystem>("structural"));
    _sys->set_eigenproblem_type(libMesh::GHEP);
    
    
    // FEType to initialize the system
    libMesh::FEType fetype (libMesh::FIRST, libMesh::LAGRANGE);
    
    // initialize the system to the right set of variables
    _sys_init    = new MAST::StructuralSystemInitialization(*_sys,
                                                               _sys->name(),
                                                               fetype);
    _discipline  = new MAST::StructuralDiscipline(*_eq_sys);
    
    // initialize the equation system
    _eq_sys->init();
    
    // set the property card data
    std::map<unsigned int, MAST::ElementPropertyCardBase*>::iterator
    it  = _elem_property_map.begin(),
    end = _elem_property_map.end();
    
    for ( ; it != end; it++)
        _discipline->set_property_for_subdomain(it->first, *it->second);
}

