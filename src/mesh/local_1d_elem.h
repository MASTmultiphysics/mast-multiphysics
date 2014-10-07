

#ifndef __mast__local_1d_elem__
#define __mast__local_1d_elem__

// MAST includes
#include "mesh/local_elem_base.h"


namespace MAST {
    
    
    /*!
     *   class provides a simple mechanism to create a geometric element
     *   in the local coordinate system.
     */
    class Local1DElem : public MAST::LocalElemBase {
    public:
        /*!
         *    constructor takes a reference to the element for this which this
         *    local element is to be created, and a vector in the element
         *    x-y plane. The x-axis, by default, is along the element length.
         */
        Local1DElem(const libMesh::Elem& elem,
                    const libMesh::Point& y);
        
        
        virtual ~Local1DElem();
        
        
        /*!
         *   Calculates the surface normal of the element at the given point
         *   \par p, where the finite element provided surface normal is
         *   \p n_local. The calculated surface normal in the global
         *   coordinate system is \p n_global. For 1D or 2D elements, the flux
         *   loads can be defined over the entire domain of the element. This,
         *   of course, is an idealization of very thin surfaces that are
         *   exchanging heat with its surroundings.
         *
         *   This method is needed because 1D and 2D elements can live in global
         *   3D space. To deal with this a geometric element with a local
         *   coordinate system is defined and used for initialization of
         *   finite elements.
         */
        virtual void domain_surface_normal_in_global_coordinates(const libMesh::Point& p,
                                                                 RealVector3& n_global) const;

    protected:
        
        /*!
         *    creation of an element in the local coordinate system
         */
        void _create_local_elem();
        
        /*!
         *    orientation of element local y-axis
         */
        libMesh::Point _local_y;
    };
    
}

#endif // __mast__local_1d_elem__
