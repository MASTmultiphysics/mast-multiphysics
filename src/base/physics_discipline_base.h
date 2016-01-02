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
    template <typename T> class ConstantFunction;
    class SystemInitialization;
    class Parameter;
    class PointLoadCondition;
    class OutputFunctionBase;
    
    
    // typedefs
    typedef std::multimap<libMesh::boundary_id_type, MAST::BoundaryConditionBase*>  SideBCMapType;
    typedef std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*> VolumeBCMapType;
    typedef std::map<libMesh::subdomain_id_type, const MAST::ElementPropertyCardBase*>      PropertyCardMapType;
    typedef std::map<libMesh::boundary_id_type, MAST::DirichletBoundaryCondition*>  DirichletBCMapType;
    typedef std::set<MAST::PointLoadCondition*> PointLoadSetType;
    typedef std::multimap<libMesh::subdomain_id_type, MAST::OutputFunctionBase*> VolumeOutputMapType;
    
    
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
         *    @returns a const reference to the volume outputs
         */
        const MAST::VolumeOutputMapType& volume_output() const{
            return _vol_output_map;
        }

        
        /*!
         *    @returns a  reference to the volume outputs
         */
        MAST::VolumeOutputMapType& volume_output() {
            
            return _vol_output_map;
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
         *    adds the output to this discipline for evaluation
         */
        void add_volume_output(libMesh::subdomain_id_type bid,
                               MAST::OutputFunctionBase& output);
        
        
        /*!
         *    initializes the system for dirichlet boundary conditions
         */
        void init_system_dirichlet_bc(libMesh::System& sys) const;

        
        /*!
         *    clears the system dirichlet boundary conditions
         */
        void clear_system_dirichlet_bc(libMesh::System& sys) const;

        
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
         *    get property card for the specified subdomain id \par i
         */
        const MAST::ElementPropertyCardBase& get_property_card(const unsigned int i) const;
        
        /*!
         *   Adds the parameter and function pairing
         */
        void add_parameter(MAST::Parameter& f);
        
        /*!
         *   Removes the parameter and function pairing
         */
        void remove_parameter(MAST::Parameter& f);
        

        /*!
         *   Returns the function corresponding to a parameter
         */
        const MAST::FunctionBase* get_parameter(const Real* par) const;
        
        
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
         *   map of sensitivity parameters and the corresponding functions that
         *   are directly dependent on these parameters
         */
        std::map<const Real*, const MAST::FunctionBase*> _parameter_map;
        
        /*!
         *   side boundary condition map of boundary id and load
         */
        MAST::SideBCMapType _side_bc_map;
        
        /*!
         *   Dirichlet boundary condition map of boundary id and load
         */
        MAST::DirichletBCMapType _dirichlet_bc_map;
        
        /*!
         *   volume boundary condition map of boundary id and load
         */
        MAST::VolumeBCMapType _vol_bc_map;
        
        /*!
         *   point loads
         */
        MAST::PointLoadSetType _point_loads;
        
        /*!
         *   volume output functions
         */
        MAST::VolumeOutputMapType _vol_output_map;
    };
    
}


#endif // __mast__physics_discipline_base__
