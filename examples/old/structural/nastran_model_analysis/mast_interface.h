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

#ifndef __mast_model_interface_h__
#define __mast_model_interface_h__

// C++ includes
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <boost/bimap.hpp>

// MAST includes
#include "base/mast_data_types.h"


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

        struct SubCase {
            SubCase(unsigned int s,
                    unsigned int l,
                    unsigned int b):
            sid  (s),
            load (l),
            spc  (b) {}

            SubCase(const MAST::Model::SubCase& sub):
            sid  (sub.sid),
            load (sub.load),
            spc  (sub.spc) {}

            unsigned int
            sid,
            load,
            spc;
        };

        struct SPC {
            SPC(unsigned int n,
                unsigned int c,
                Real         v):
            node    (n),
            comp    (c),
            val     (v) {}
            
            SPC(const MAST::Model::SPC& spc):
            node  (spc.node),
            comp  (spc.comp),
            val   (spc.val) {}
            
            unsigned int
            node,
            comp,
            val;
        };

        struct Force {
            Force(unsigned int n,
                  Real         fx,
                  Real         fy,
                  Real         fz):
            node    (n),
            f_x     (fx),
            f_y     (fy),
            f_z     (fz)    {}
            
            Force(const MAST::Model::Force& force):
            node  (force.node),
            f_x   (force.f_x),
            f_y   (force.f_y),
            f_z   (force.f_z) {}
            
            unsigned int
            node,
            f_x,
            f_y,
            f_z;
        };

        
        Model(libMesh::LibMeshInit& init, const std::string& nm);
        
        ~Model();
        
        void
        set_sol(unsigned int sol) {_sol = sol;}
        
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

        void
        add_subcase(int sid, int load_set, int spc);

        void
        add_spc(int node, int component, Real val);

        void
        add_force(int node, Real fx, Real fy, Real fz);

        
        
        const std::string&           file_path()       {return _file;}
        
        unsigned int                 get_sol()         {return _sol;}
        
        unsigned int                 n_subcases()      {return (unsigned int)_subcases.size();}

        const std::vector<MAST::Model::SubCase>&  get_subcases() const {return _subcases;}
        
        const MAST::Model::SubCase&  get_subcase(unsigned int i) const;
        
        libMesh::MeshBase&           get_mesh()        { return *_mesh; }
        
        libMesh::EquationSystems&    get_eq_sys()      { return *_eq_sys; }

        MAST::NonlinearSystem&       get_system()      { return *_sys; }

        MAST::SystemInitialization&  get_system_init() { return *_sys_init; }
        
        MAST::PhysicsDisciplineBase& get_discipline()  { return *_discipline; }
        
        std::vector<MAST::Model::SPC>& get_spcs()      { return _spcs;}
        
        std::vector<MAST::Model::Force>& get_forces()  { return _forces;}
        
        boost::bimap<unsigned int, libMesh::Node*>&
        get_node_id_map() {return _node_id_map;}

        void
        initialize_after_mesh();
        
        
        void
        clear_loads();
        
        void
        print(std::ostream& o);
        
    protected:
        
        /*!
         *   name of file
         */
        std::string         _file;
        
        
        /*!
         *  mesh data structure
         */
        libMesh::SerialMesh  *_mesh;
        
        
        /*!
         *  solution type in the file
         */
        unsigned int            _sol;
        
        
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
        boost::bimap<unsigned int, libMesh::Node*> _node_id_map;
        
        /*!
         *  map of elem IDs to pointers
         */
        boost::bimap<unsigned int, libMesh::Elem*> _elem_id_map;

        
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

        /*!
         *  vector of subcases to be solved
         */
        std::vector<MAST::Model::SubCase>  _subcases;

        /*!
         *  vector of subcases to be solved
         */
        std::vector<MAST::Model::SPC>     _spcs;

        /*!
         *  vector of subcases to be solved
         */
        std::vector<MAST::Model::Force>   _forces;
    };
}


#endif // __mast_model_interface_h__

