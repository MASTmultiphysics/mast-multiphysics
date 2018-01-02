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

#ifndef __mast__flutter_root_base_h__
#define __mast__flutter_root_base_h__

// MAST includes
#include "base/mast_data_types.h"


namespace MAST {

    class FlutterRootBase {
        
    public:
        
        /*!
         *  default constructor
         */
        FlutterRootBase();
        
        /*!
         *   copy constructor
         */
        FlutterRootBase(const MAST::FlutterRootBase& f);
        
        virtual void copy_root(const MAST::FlutterRootBase& f);
        
        virtual ~FlutterRootBase() {}
        
        bool has_sensitivity_data, if_nonphysical_root;
        
        Real kr, g, kr_sens, V, omega, V_sens;
        
        Complex root, root_sens;
        
        /*!
         *    right and left eigenvevtors
         */
        ComplexVectorX  eig_vec_right, eig_vec_left;
        
        RealVectorX modal_participation;
        
    };

}


#endif // __mast__flutter_root_base_h__
