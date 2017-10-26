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

#ifndef __mast_nastran_model_analysis_h__
#define __mast_nastran_model_analysis_h__



// MAST includes
#include "base/mast_data_types.h"


// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/fe_type.h"
#include "libmesh/dof_map.h"


namespace MAST {
    
    // Forward declerations
    class Model;
    class Parameter;
    class StructuralSystemInitialization;
    class StructuralDiscipline;
    class DirichletBoundaryCondition;
    class BoundaryConditionBase;
    class StressStrainOutputBase;
    class NonlinearSystem;

    
    class NastranModelAnalysis {
        
    public:
        
        NastranModelAnalysis();
        
        
        ~NastranModelAnalysis();
        
                
        /*!
         *   initializes the object for specified characteristics
         */
        void init(libMesh::ElemType etype, bool if_nonlin);
        

        MAST::Parameter* get_parameter(const std::string& nm);

        
        /*!
         *  solves the system and returns the final solution
         */
        const libMesh::NumericVector<Real>&
        solve(bool if_write_output = false);
        
        
        /*!
         *  solves the sensitivity of system and returns the final solution
         */
        const libMesh::NumericVector<Real>&
        sensitivity_solve(MAST::Parameter& p,
                          bool if_write_output = false);

        
    protected:

        /*!
         *   solves the current subcase for which BC and loads have been
         *   imported into the model
         */
        void _solve_subcase();

        /*!
         *   initializes the dof constraints
         */
        void _initialize_dof_constraints();

        
        /*!
         *   linear static analysis
         */
        void _linear_static_solve();

        /*!
         *   modal analysis
         */
        void _normal_modes_solve();
        
        /*!
         *   FEM model
         */
        MAST::Model* _model;
        
    };
}


#endif //__mast_nastran_model_analysis_h__
