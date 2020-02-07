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

#ifndef __mast__gaf_database_h__
#define __mast__gaf_database_h__

// C++ includes
#include <vector>
#include <map>

// MAST includes
#include "base/mast_data_types.h"
#include "elasticity/fsi_generalized_aero_force_assembly.h"

// libMesh includes
#include "libmesh/numeric_vector.h"


namespace MAST{
    
    // Forward decleraitons
    class Parameter;
    class FrequencyFunction;
    
    class GAFDatabase:
    public MAST::FSIGeneralizedAeroForceAssembly {
        
    public:
        
        GAFDatabase(const unsigned int n_modes);
        
        ~GAFDatabase() { }
        
        void init(MAST::FrequencyFunction*                        freq,
                  MAST::ComplexSolverBase*                complex_solver,
                  MAST::PressureFunction*                 pressure_func,
                  MAST::FrequencyDomainPressureFunction*  freq_pressure_func,
                  MAST::ComplexMeshFieldFunction*         displ_func);

        
        void
        set_evaluate_mode(bool f);
        
        
        void
        write_gaf_file(const std::string& nm,
                       std::vector<libMesh::NumericVector<Real>*>& modes);
        
        
        
        void
        read_gaf_file(const std::string& nm,
                      std::vector<libMesh::NumericVector<Real>*>& modes);
        
        
        ComplexMatrixX&
        add_kr_mat(const Real kr,
                   const ComplexMatrixX& mat,
                   const bool if_kr_sens);
        
        ComplexMatrixX
        get_kr_mat(const Real kr,
                   const std::map<Real, ComplexMatrixX>& data);
        
        
        virtual void
        assemble_generalized_aerodynamic_force_matrix
        (std::vector<libMesh::NumericVector<Real>*>& basis,
         ComplexMatrixX& mat,
         MAST::Parameter* p = nullptr);
        
        
        
    protected:

        MAST::FrequencyFunction*                    _freq;
        bool                                _if_evaluate;
        unsigned int                        _n_modes;
        std::map<Real, ComplexMatrixX>      _kr_to_gaf_map;
        std::map<Real, ComplexMatrixX>      _kr_to_gaf_kr_sens_map;
    };
}


#endif  // __mast__gaf_database_h__

