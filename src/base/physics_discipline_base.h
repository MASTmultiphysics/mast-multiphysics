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

#ifndef __mast__physics_discipline_base__
#define __mast__physics_discipline_base__

// C++ includes
#include <map>

// MAST includes
#include "base/mast_data_types.h"

// libMesh includes
#include "libmesh/equation_systems.h"


namespace MAST {
    
    // Forward declerations
    class BoundaryConditionBase;
    class DirichletBoundaryCondition;
    class PropertyCardBase;
    class FunctionBase;
    class FunctionSetBase;
    class ElementPropertyCardBase;
    class SystemInitialization;
    class Parameter;
    class PointLoadCondition;
    class NonlinearSystem;
    class GeomElem;
    
    
    // typedefs
    typedef std::multimap<libMesh::boundary_id_type, MAST::BoundaryConditionBase*>  SideBCMapType;
    typedef std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*> VolumeBCMapType;
    typedef std::map<libMesh::subdomain_id_type, const MAST::ElementPropertyCardBase*>      PropertyCardMapType;
    typedef std::map<libMesh::boundary_id_type, MAST::DirichletBoundaryCondition*>  DirichletBCMapType;
    typedef std::set<MAST::PointLoadCondition*> PointLoadSetType;
    
    class PhysicsDisciplineBase {
    public:
        
        
        // Constructor
        PhysicsDisciplineBase(libMesh::EquationSystems& eq_sys):
        _eq_systems(eq_sys)
        { }
        
        /*!
         *   virtual destructor
         */
        virtual ~PhysicsDisciplineBase()
        { }
        
        
        /*!
         *   returns a reference to the libMesh::System object
         */
        libMesh::EquationSystems& get_equation_systems() {
            return _eq_systems;
        }
        
        /*!
         *   clear the loads and pointer to static solution system for
         *   this structural model
         */
        void clear_loads();
        
        /*!
         *   clear the specified volume load from the applied loads
         */
        void clear_volume_load(libMesh::subdomain_id_type sid,
                               MAST::BoundaryConditionBase& load);
        
        /*!
         *   adds the specified side loads for the boudnary with tag \p b_id
         */
        void add_side_load(libMesh::boundary_id_type bid,
                           MAST::BoundaryConditionBase& load);
        
        /*!
         *  adds the specified Dirichlet boundary condition for the boundary
         *  with tag \p b_id
         */
        void add_dirichlet_bc(libMesh::boundary_id_type bid,
                              MAST::DirichletBoundaryCondition& load);
        
        /*!
         *    @returns a const reference to the side boundary conditions
         */
        const MAST::SideBCMapType& side_loads() const {
            return _side_bc_map;
        }
        
        
        /*!
         *    @returns a reference to the side boundary conditions
         */
        MAST::SideBCMapType& side_loads() {
            return _side_bc_map;
        }
        
        /*!
         *   adds the specified volume loads for the elements with
         *   subdomain tag \p s_id
         */
        void add_volume_load(libMesh::subdomain_id_type bid,
                             MAST::BoundaryConditionBase& load);

        
        /*!
         *   adds the specified point load
         */
        void add_point_load(MAST::PointLoadCondition& load);

        
        /*!
         *    @returns a const reference to the volume boundary conditions
         */
        const MAST::VolumeBCMapType& volume_loads() const{
            return _vol_bc_map;
        }
        
        /*!
         *    @returns a reference to the volume boundary conditions
         */
        MAST::VolumeBCMapType& volume_loads() {
            return _vol_bc_map;
        }

        
        /*!
         *    @returns a const reference to the point load boundary conditions
         */
        const MAST::PointLoadSetType& point_loads() const{
            return _point_loads;
        }
        
        /*!
         *    @returns a reference to the point load boundary conditions
         */
        MAST::PointLoadSetType& point_loads() {
            return _point_loads;
        }

        
        
        /*!
         *    constrain dofs on a subdomain to zero
         */
        void constrain_subdomain_dofs_for_var(const libMesh::subdomain_id_type sid,
                                              const unsigned int var);
        
        
        /*!
         *    initializes the system for dirichlet boundary conditions
         */
        void init_system_dirichlet_bc(MAST::NonlinearSystem& sys) const;

        
        /*!
         *    clears the system dirichlet boundary conditions
         */
        void clear_system_dirichlet_bc(MAST::NonlinearSystem& sys) const;

        
        /*!
         *   Prepares a list of the constrained dofs for system \p sys and 
         *   returns in \p dof_ids.
         */
        void get_system_dirichlet_bc_dofs(libMesh::System& sys,
                                          std::set<unsigned int>& dof_ids) const;
        
        /*!
         *    sets the same property for all elements in the specified subdomain
         */
        void set_property_for_subdomain(const libMesh::subdomain_id_type sid,
                                        const MAST::ElementPropertyCardBase& prop);
        
        /*!
         *    get property card for the specified element
         */
        const MAST::ElementPropertyCardBase& get_property_card(const libMesh::Elem& elem) const;

        /*!
         *    get property card for the specified element
         */
        const MAST::ElementPropertyCardBase& get_property_card(const MAST::GeomElem& elem) const;
        
        /*!
         *    get property card for the specified subdomain id \p i
         */
        const MAST::ElementPropertyCardBase& get_property_card(const unsigned int sid) const;
        
        
        
    protected:
        
        /*!
         *    libMesh::System for which analysis is to be performed
         */
        libMesh::EquationSystems& _eq_systems;
        
        /*!
         *   map of element property cards for each element
         */
        MAST::PropertyCardMapType _element_property;
        
        /*!
         *   side boundary condition map of boundary id and load
         */
        MAST::SideBCMapType _side_bc_map;
        
        /*!
         *   Dirichlet boundary condition map of boundary id and load
         */
        MAST::DirichletBCMapType _dirichlet_bc_map;
        
        
        /*!
         *   variables constrained on subdomain
         */
        std::map<libMesh::subdomain_id_type, std::vector<unsigned int> > _subdomain_var_constraint;
        
        /*!
         *   volume boundary condition map of boundary id and load
         */
        MAST::VolumeBCMapType _vol_bc_map;
        
        /*!
         *   point loads
         */
        MAST::PointLoadSetType _point_loads;
    };
    
}


#endif // __mast__physics_discipline_base__
