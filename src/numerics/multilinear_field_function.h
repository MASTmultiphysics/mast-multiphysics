
#ifndef __mast__multilinear_field_function__
#define __mast__multilinear_field_function__

// C++ includes
#include <map>

// MAST includes
#include "base/field_function_base.h"


namespace MAST {

    /*!
     *   the following is used for calculation of the return value
     *   f(x) is defined for x for each x0 < x < x1
     *   if   x <= x0,      f(x) = f(x0)
     *   if   x0 < x < x1,  f(x) is interpolated
     *   if   x >= x1,      f(x) = f(x1)
     */
    template <typename ValType>
    class MultilinearInterpolation:
    public MAST::FieldFunction<ValType> {
        
    public:
        MultilinearInterpolation(const std::string& nm,
                                 MAST::FieldFunction<ValType>* abscissa,
                                 std::map<Real, MAST::FieldFunction<ValType>*>& values);
        
        
        MultilinearInterpolation(const MAST::MultilinearInterpolation<ValType>& o);
        
        
        virtual ~MultilinearInterpolation();

        virtual std::auto_ptr<MAST::FieldFunction<ValType> > clone() const;
        
        
        virtual void operator() (const libMesh::Point& p, Real t, Real& v) const;
      
        
        /*!
         *   calculates derivative of the function
         */
        virtual void derivative(const MAST::DerivativeType d,
                                const MAST::FunctionBase& f,
                                const libMesh::Point& p,
                                const Real t,
                                ValType& v) const;
        
    protected:
        
        MAST::FieldFunction<ValType>* _abscissa;

        std::map<Real, MAST::FieldFunction<ValType>*> _values;
        
    };
}


#endif //__mast__multilinear_field_function__
