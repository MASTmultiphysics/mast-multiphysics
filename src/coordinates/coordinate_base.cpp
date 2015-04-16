
// MAST includes
#include "coordinates/coordinate_base.h"


MAST::CoordinateBase::CoordinateBase(const std::string& nm):
MAST::FieldFunction<RealMatrixX>(nm) {
    
}



void
MAST::CoordinateBase::
stress_strain_transformation_matrix(const RealMatrixX& T,
                                    RealMatrixX &mat) {
    
    mat.setZero(6,6);
    
    std::vector<std::pair<unsigned int, unsigned int> > ids(6);
    ids[0] = std::pair<unsigned int, unsigned int>(0,0);
    ids[1] = std::pair<unsigned int, unsigned int>(1,1);
    ids[2] = std::pair<unsigned int, unsigned int>(2,2);
    ids[3] = std::pair<unsigned int, unsigned int>(0,1);
    ids[4] = std::pair<unsigned int, unsigned int>(1,2);
    ids[5] = std::pair<unsigned int, unsigned int>(2,0);
    
    for (unsigned int i=0; i<6; i++) {
        for (unsigned int j=0; j<3; j++) {
            // the first three columns of the matrix
            mat(i,j)   =
            T(ids[j].first, ids[i].first) * T(ids[j].second, ids[i].second);
            //            std::cout
            //            << "( " << i << ", " << j << " ): "
            //            << ids[j].first+1 << ids[i].first+1 << " x "
            //            << ids[j].second+1 << ids[i].second+1 << std::endl;
            
            // last three columns of the matrix
            mat(i,j+3) =
            T(ids[j+3].first, ids[i].first)  * T(ids[j+3].second, ids[i].second) +
            T(ids[j+3].first, ids[i].second) * T(ids[j+3].second, ids[i].first);
            
            //            std::cout
            //            << "( " << i << ", " << j+3 << " ): "
            //            << ids[j+3].first+1 << ids[i].first+1  << " x "
            //            <<  ids[j+3].second+1 << ids[i].second+1
            //            << " + " << ids[j+3].first+1 << ids[i].second+1 << " x "
            //            << ids[j+3].second+1 << ids[i].first+1 << std::endl;
        }
        std::cout << std::endl;
    }
}

