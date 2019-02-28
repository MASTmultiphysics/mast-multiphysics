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

#ifndef __mast__constrain_beam_dofs__
#define __mast__constrain_beam_dofs__

// C++ includes
#include <vector>


// libMesh includes
#include "libmesh/system.h"

namespace MAST {
    
    // Forward declerations
    class NonlinearSystem;
    
    
    class ConstrainBeamDofs:
    public libMesh::System::Constraint  {
        
    public:
        ConstrainBeamDofs(MAST::NonlinearSystem& sys);
        
        virtual ~ConstrainBeamDofs();
        
        /*!
         *   This will constrain all dofs for a variable. This is helpful for
         *   eigenvalue anlaysis when one set of deformation is unnecessary and
         *    its presence can complicate the analysis.
         */
        virtual void constrain();
        
    protected:

        MAST::NonlinearSystem& _system;
    };
}

#endif // __mast__constrain_beam_dofs__

