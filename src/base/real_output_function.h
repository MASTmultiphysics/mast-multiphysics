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

#ifndef __mast__real_output_function_h__
#define __mast__real_output_function_h__

// C++ includes
#include <map>

// MAST includes
#include "base/mast_data_types.h"
#include "base/output_function_base.h"
#include "base/physics_discipline_base.h"



namespace MAST {
    
    
    class RealOutputFunction:
    public  MAST::OutputFunctionBase {
        
    public:
        
        /*!
         *   default constructor
         */
        RealOutputFunction(MAST::OutputQuantityType t);
        
        
        /*!
         *   destructir
         */
        virtual ~RealOutputFunction();
        
        
        /*!
         *   zeros all data and removes sensitivity or derivative information
         */
        void clear();
        
        
        void set_value(Real val);


        void add_value(Real val);

        
        Real get_value() const;
        
        
        void set_derivative(const RealVectorX& dvdX);
        
        
        const RealVectorX& get_derivative() const;
        
        
        void add_sensitivity(const MAST::FunctionBase* f,
                             const Real dvdp);

        
        void set_sensitivity(const MAST::FunctionBase* f,
                             const Real dvdp);
        
        
        Real get_sensitivity(const MAST::FunctionBase* f) const;
        
        
    protected:
        
        Real _value;

        RealVectorX _derivative;
        
        std::map<const MAST::FunctionBase*, Real> _sensitivity;
    };
}



#endif // __mast__real_output_function_h__
