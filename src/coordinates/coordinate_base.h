

#ifndef __mast__coordinate_base__
#define __mast__coordinate_base__


// MAST includes
#include "base/field_function_base.h"


namespace MAST {
    
    /*!
     *    Provides the transformation matrix T to transform
     *    vector from the orientation provided in this matrix,
     *    to one in the global basis
     */
    class CoordinateBase:
    public MAST::FieldFunction<RealMatrixX> {
        
    public:
        
        CoordinateBase(const std::string& nm);

        /*!
         *   prepares the matrix \mat that transforms stress and strain tensors
         *   represented in a 6x1 vector from the coordinate system in _orient
         *   to the global coordinate system. Note that the shear straints in
         *   the strain tensor vector should be represented in the tensor quantities,
         *   and not the engineering strain.
         */
        void stress_strain_transformation_matrix(const RealMatrixX& T,
                                                 RealMatrixX& mat);
        

    protected:
        
    };
}




#endif // __mast__coordinate_base__