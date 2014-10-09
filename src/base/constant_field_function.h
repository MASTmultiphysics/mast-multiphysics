
#ifndef __mast__constant_field_function__
#define __mast__constant_field_function__

// MAST includes
#include "base/field_function_base.h"


namespace MAST {
    
    // Forward declerations
    class Parameter;
    
    class ConstantFieldFunction:
    public MAST::FieldFunction<Real> {
    
    public:
        
        ConstantFieldFunction(const std::string& nm,
                              const MAST::Parameter& p);

        
        ConstantFieldFunction(const MAST::ConstantFieldFunction& f);

        
        virtual ~ConstantFieldFunction();

        /*!
         *   @returns a clone of the function
         */
        virtual std::auto_ptr<MAST::FieldFunction<Real> > clone() const;
        
        
        /*!
         *    calculates the value of the function at the specified point,
         *    \par p, and time, \par t, and returns it in \p v.
         */
        virtual void operator() (const libMesh::Point& p,
                                 const Real t,
                                 Real& v) const;
        
        
        /*!
         *    calculates the value of the function at the specified point,
         *    \par p, and time, \par t, and returns it in \p v.
         */
        virtual void derivative (const MAST::DerivativeType d,
                                 const MAST::FunctionBase& f,
                                 const libMesh::Point& p,
                                 const Real t,
                                 Real& v) const;

        
        
    protected:

        /*!
         *   parameter which defines this field function
         */
        const MAST::Parameter& _p;
        
    };
}


#endif // __mast__constant_field_function__
