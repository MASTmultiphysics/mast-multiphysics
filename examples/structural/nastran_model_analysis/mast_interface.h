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

#ifndef __mast_model_interface_h__
#define __mast_model_interface_h__

// C++ includes
#include <map>
#include <string>


// libMesh includes
#include "libmesh/equation_systems.h"
#include "libmesh/serial_mesh.h"
#include "libmesh/node.h"
#include "libmesh/elem.h"


namespace MAST {
    
    
    // Forward decleration
    class Parameter;
    class FunctionBase;
    class MaterialPropertyCardBase;
    class ElementPropertyCardBase;
    class NonlinearSystem;
    class SystemInitialization;
    class PhysicsDisciplineBase;
    

    /*!
     *    This class defines an interface to conveniently create and store 
     *    models. The motivation for this class came from importing pyNastran
     *    models. Some assumptions are made for simplicity here:
     *
     *     - PSHELL allows for different materials for membrane and bending. This
     *         class assumes the same material for membrane and bending.
     *     - the shear correction factor is automatically added to material property
     */
    class Model {
        
    public:
        
        Model(libMesh::LibMeshInit& init);
        
        ~Model();
        
        void
        add_node(int n_id, double x, double y, double z);
        
        void
        add_edge2(int e_id, int p_id,   int n1,   int n2);

        void
        add_quad4(int e_id, int p_id,   int n1,   int n2,   int n3,   int n4);

        void
        add_isotropic_material(int m_id, double E, double nu, double rho, double alpha);

        void
        add_2d_section_property(int p_id, int m_id, double t);
        
        libMesh::MeshBase&           get_mesh()        { return *_mesh; }
        
        libMesh::EquationSystems&    get_eq_sys()      { return *_eq_sys; }

        MAST::NonlinearSystem&       get_system()      { return *_sys; }

        MAST::SystemInitialization&  get_system_init() { return *_sys_init; }
        
        MAST::PhysicsDisciplineBase& get_discipline()  { return *_discipline; }

        void
        initialize_after_mesh();
        
    protected:
        
        /*!
         *  mesh data structure
         */
        libMesh::SerialMesh  *_mesh;
        
        /*!
         *   equation systems
         */
        libMesh::EquationSystems *_eq_sys;
        
        /*!
         *   system
         */
        MAST::NonlinearSystem *_sys;

        /*!
         *   system initializer
         */
        MAST::SystemInitialization *_sys_init;
        
        /*!
         *   discipline
         */
        MAST::PhysicsDisciplineBase *_discipline;
        
        /*!
         *  map of node IDs to pointers
         */
        std::map<unsigned int, libMesh::Node*> _node_id_map;
        
        /*!
         *  map of elem IDs to pointers
         */
        std::map<unsigned int, libMesh::Elem*> _elem_id_map;

        
        /*!
         *  map of parameters
         */
        std::map<std::string, MAST::Parameter*> _param_map;

        /*!
         *  map of parameters
         */
        std::map<std::string, MAST::FunctionBase*> _function_map;

        /*!
         *  map of material properties
         */
        std::map<unsigned int, MAST::MaterialPropertyCardBase*> _material_map;

        /*!
         *  map of element section properties
         */
        std::map<unsigned int, MAST::ElementPropertyCardBase*> _elem_property_map;

        
    };
}


#endif // __mast_model_interface_h__

