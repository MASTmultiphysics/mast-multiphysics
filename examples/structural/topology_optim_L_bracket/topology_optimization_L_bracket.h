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

#ifndef __mast_topology_optimization_level_set_L_bracket_h__
#define __mast_topology_optimization_level_set_L_bracket_h__


// MAST includes
#include "examples/structural/topology_optim_2D/topology_optim_2D.h"


namespace MAST  {
    
    
    namespace Examples {
        
        class TopologyOptimizationLevelSetLBracket:
        public MAST::Examples::TopologyOptimizationLevelSet2D {
            
        public:
            
            TopologyOptimizationLevelSetLBracket(const libMesh::Parallel::Communicator& comm_in);
            
            virtual ~TopologyOptimizationLevelSetLBracket();
            
            
        protected:
            
            virtual void _init_mesh();
            virtual void _init_dirichlet_conditions();
            virtual void _init_loads();
            virtual void _init_phi_dvs();
            virtual void _delete_elems_from_mesh(libMesh::MeshBase& mesh);
        };
    }
}


#endif //  __mast_topology_optimization_level_set_L_bracket_h__

