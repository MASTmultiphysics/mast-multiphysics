

#ifndef __mast__basis_matrix_coordinate__
#define __mast__basis_matrix_coordinate__


// MAST includes
#include "coordinates/coordinate_base.h"


namespace MAST {
    
    /*!
     *    Provides the transformation matrix T to transform
     *    vector from the orientation provided in this matrix,
     *    to one in the global basis
     */
    class BasisMatrixCoordinate:
    public MAST::FieldFunction<RealMatrixX> {
        
    public:
        BasisMatrixCoordinate(const std::string& nm,
                              MAST::FieldFunction<RealMatrixX>* basis);
        
        virtual ~BasisMatrixCoordinate();
        
        /*!
         *  Copy contructor
         */
        BasisMatrixCoordinate(const MAST::BasisMatrixCoordinate& c);
        
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
        
        MAST::FieldFunction<RealMatrixX>* _basis;
        
    };
}




#endif // __mast__basis_matrix_coordinate__

