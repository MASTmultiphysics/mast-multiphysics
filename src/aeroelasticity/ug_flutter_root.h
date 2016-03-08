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

#ifndef __mast__ug_flutter_root_h__
#define __mast__ug_flutter_root_h__

// MAST includes
#include "aeroelasticity/flutter_root_base.h"


namespace MAST {
    
    class UGFlutterRoot:
    public MAST::FlutterRootBase {
    
    public:
        
        /*!
         *  default constructor
         */
        UGFlutterRoot();
        
        /*!
         *   copy constructor
         */
        UGFlutterRoot(const MAST::UGFlutterRoot& f);
        
        
        virtual ~UGFlutterRoot() {}

        
        /*!
         *   initializes the data
         */
        virtual void init(const Real kr_ref_val,
                          const Real b_ref,
                          const Complex num,
                          const Complex den,
                          const ComplexMatrixX& Bmat,
                          const ComplexVectorX& evec_right,
                          const ComplexVectorX& evec_left);
        
        
    };
}


#endif // __mast__ug_flutter_root_h__
