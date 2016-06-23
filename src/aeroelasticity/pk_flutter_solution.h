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

#ifndef __mast__pk_flutter_solution_h__
#define __mast__pk_flutter_solution_h__


// C++ includes
#include <vector>


// MAST includes
#include "aeroelasticity/flutter_solution_base.h"


namespace MAST {
    
    // Forward declerations
    class PKFlutterSolver;
    class LAPACK_ZGGEV;
    
    
    class PKFlutterSolution:
    public MAST::FlutterSolutionBase {
    public:
        
        PKFlutterSolution():
        MAST::FlutterSolutionBase()
        { }
        
        
        virtual ~PKFlutterSolution() {}
        
        /*!
         *   initializes the flutter solution from an eigensolution
         */
        virtual void init (const MAST::PKFlutterSolver& solver,
                           const Real k_red,
                           const Real v_ref,
                           const Real bref,
                           const RealMatrixX& kmat,
                           const MAST::LAPACK_ZGGEV& eig_sol);

    protected:
        
        /*!
         *  value of reduced frequency for this solution
         */
        Real  _k_red;
        
        /*!
         *   structural stiffness matrix
         */
        RealMatrixX _stiff_mat;

    };
}

#endif // __mast__pk_flutter_solution_h__

