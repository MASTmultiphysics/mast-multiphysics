

#ifndef __mast__local_3d_elem__
#define __mast__local_3d_elem__


// MAST includes
#include "mesh/local_elem_base.h"


namespace MAST {
    
    /*!
     *   class provides a simple mechanism to create a geometric element
     *   in the local coordinate system.
     */
    class Local3DElem:
    public MAST::LocalElemBase {
        
    public:
        Local3DElem(const libMesh::Elem& elem);
        
        
        virtual ~Local3DElem();
        
        
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
        virtual void
        domain_surface_normal_in_global_coordinates(const libMesh::Point& p,
                                                    RealVector3& n_global) const;
        
    protected:
        
    };
    
}


#endif // __mast__local_3d_elem__
