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

#ifndef __mast__time_domain_flutter_root_base_h__
#define __mast__time_domain_flutter_root_base_h__

// MAST includes
#include "base/mast_data_types.h"


namespace MAST {

    class TimeDomainFlutterRootBase {
        
    public:
        
        /*!
         *  default constructor
         */
        TimeDomainFlutterRootBase();
        
        /*!
         *   copy constructor
         */
        TimeDomainFlutterRootBase(const TimeDomainFlutterRootBase& f);


        /*!
         *   initializes the data
         */
        void init(const Real v_ref_val,
                  const Complex num,
                  const Complex den,
                  const RealMatrixX& Bmat,
                  const ComplexVectorX& evec_right,
                  const ComplexVectorX& evec_left);

        
        void copy_root(const MAST::TimeDomainFlutterRootBase& f);
        
        virtual ~TimeDomainFlutterRootBase() {}
        
        bool has_sensitivity_data;
        
        Real V, omega, V_sens;
        
        Complex root, root_sens;
        
        /*!
         *    right and left eigenvevtors
         */
        ComplexVectorX  eig_vec_right, eig_vec_left;
        
        RealVectorX modal_participation;
        
    };

}


#endif // __mast__time_domain_flutter_root_base_h__
