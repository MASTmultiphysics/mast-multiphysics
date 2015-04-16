

#ifndef __mast__polar_coordinate_base__
#define __mast__polar_coordinate_base__


// MAST includes
#include "coordinates/coordinate_base.h"


namespace MAST {
    
    /*!
     *    Defines a polar coordinate system with the radius and 
     *    obtained from the two parameters provided in the constructor
     */
    class PolarCoordinate:
    public MAST::CoordinateBase {
        
    public:
        PolarCoordinate(const std::string& nm,
                        MAST::FieldFunction<Real>* theta);
        
        virtual ~PolarCoordinate();
        
        /*!
         *  Copy contructor
         */
        PolarCoordinate(const MAST::PolarCoordinate& c);
        
        /*!
         *   @returns a clone of the function
         */
        virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> > clone() const;
        
        /*!
         *    calculates the value of the function at the specified point,
         *    \par p, and time, \par t, and returns it in \p v.
         */
        virtual void operator() (const libMesh::Point& p,
                                 const Real t,
                                 RealMatrixX& v) const;
        
        
        /*!
         *    calculates the value of the function at the specified point,
         *    \par p, and time, \par t, and returns it in \p v.
         */
        virtual void derivative (const MAST::DerivativeType d,
                                 const MAST::FunctionBase& f,
                                 const libMesh::Point& p,
                                 const Real t,
                                 RealMatrixX& v) const;
        
        
    protected:
        
        MAST::FieldFunction<Real>* _theta;
    };
}




#endif // __mast__polar_coordinate_base__
