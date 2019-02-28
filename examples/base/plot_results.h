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

#ifndef __mast__plot_results_h__
#define __mast__plot_results_h__

// C++ includes
#include <string>

// MAST includes
#include "base/mast_data_types.h"

// libMesh includes
#include "libmesh/system.h"



namespace MAST {
    
    // forward declerations
    class ComplexMeshFieldFunction;
    class ComplexSolverBase;

    void
    plot_structural_flutter_solution(const std::string& nm,
                                     libMesh::System& sys,
                                     const ComplexVectorX& eig_vec,
                                     const std::vector<libMesh::NumericVector<Real>*>& basis);

    void
    plot_fluid_flutter_solution(const std::string&                      nm,
                                libMesh::System&            structural_sys,
                                libMesh::System&                 fluid_sys,
                                MAST::ComplexMeshFieldFunction&      displ,
                                MAST::ComplexSolverBase&            solver,
                                const ComplexVectorX&               eig_vec,
                                const std::vector<libMesh::NumericVector<Real>*>& structural_basis);

}



#endif // __mast__plot_results_h__
