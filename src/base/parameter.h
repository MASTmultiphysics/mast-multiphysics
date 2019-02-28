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

#ifndef __mast__parameter__
#define __mast__parameter__


// MAST includes
#include "base/function_base.h"


namespace MAST {
    
    /*!
     *    This is a function that does not change.
     */
    class Parameter:
    public MAST::FunctionBase {

    public:
        
        Parameter(const std::string& nm,
                  const Real& val):
        MAST::FunctionBase(nm, false),
        _val(new Real)
        { *_val = val;}

        
        Parameter(const MAST::Parameter& f):
        MAST::FunctionBase(f),
        _val(f._val)
        { }
        

        ~Parameter() {

            delete _val;
        }
        
        
        /*!
         *   @returns a writable reference to this parameter value
         */
        Real& operator() () {
            return *_val;
        }

        
        /*!
         *   @returns the parameter value
         */
        Real operator() () const {
            return *_val;
        }
        
        
        /*!
         *    @returns the pointer to value of this function.
         */
        Real* ptr()
        {   return _val; }
        
        
        
        /*!
         *  @returns  \p true only if the given function is this.
         */
        virtual bool depends_on(const MAST::FunctionBase& f) const {
            if (&f == this)
                return true;
            else
                return false;
        }
        
        
        
        /*!
         *  sets the value of this function
         */
        void operator =(const Real& val)
        {   *_val = val; }
        
        
    protected:
        
        /*!
         *    Pointer to the value of the parameter
         */
        Real* _val;
    };
}

#endif // __mast__parameter__
