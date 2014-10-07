
#ifndef __mast__local_elem_base__
#define __mast__local_elem_base__


// MAST includes
#include "base/mast_data_types.h"


// libMesh includes
#include "libmesh/elem.h"



namespace MAST {

    // Forward declerations
    template <typename ValType> class FieldFunction;
    
    /*!
     *   class provides a simple mechanism to create a geometric element
     *   in the local coordinate system.
     */
    class LocalElemBase {
    public:
        LocalElemBase(const libMesh::Elem& elem):
        _elem(elem),
        _local_elem(NULL)
        { }
        
        
        virtual ~LocalElemBase();
        
        /*!
         *   returns a constant reference to the global element.
         */
        const libMesh::Elem& global_elem() const {
            return _elem;
        }
        
        
        /*!
         *   returns a constant reference to the local element.
         */
        const libMesh::Elem& local_elem() const {
            if (!_local_elem) // original element lies in the xy-plane
                return _elem;
            else
                return *_local_elem;
        }
        
        /*!
         *    returns the transformation matrix for this element. This is used
         *    to map the coordinates from local to global coordinate system
         */
        const RealMatrix3& T_matrix() const {
            return _T_mat;
        }
        
        
        /*!
         *    returns the transformation matrix for this element. This is used
         *    to map the coordinates from local to global coordinate system
         */
        std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
        T_matrix_function() const;
        
        
        /*!
         *    Local elements are defined for 1D and 2D elements that exist in
         *    3D space. These elements have a local coordinate system associated
         *    with the local coordinate. This method accepts the point defined
         *    the local coordinate system as the input and maps it to the
         *    global coordinate system.
         */
        void global_coordinates_location(const libMesh::Point& local,
                                         libMesh::Point& global) const;

        
        /*!
         *   Calculates the side surface normal of the element at the given
         *   point \par p, where the finite element provided surface normal is
         *   \p n_local. The calculated surface normal in the global
         *   coordinate system is \p n_global.
         *
         *   This method is needed because 1D and 2D elements can live in global
         *   3D space. To deal with this a geometric element with a local
         *   coordinate system is defined and used for initialization of
         *   finite elements.
         */
        void global_coordinates_normal(const libMesh::Point& local,
                                       RealVector3& global) const;


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
                                                                 RealVector3& n_global) const = 0;

    protected:
        
        /*!
         *   given element in global coordinate system
         */
        const libMesh::Elem& _elem;
        
        /*!
         *   element created in local coordinate system
         */
        libMesh::Elem* _local_elem;
        
        /*!
         *   nodes for local element
         */
        std::vector<libMesh::Node*> _local_nodes;
        
        /*!
         *    Transformation matrix defines T_ij = V_i^t . Vn_j, where
         *    V_i are the unit vectors of the global cs, and Vn_j are the
         *    unit vectors of the local cs. To transform a vector from global to
         *    local cs,    an_j = T^t a_i, and the reverse transformation is
         *    obtained as  a_j  = T  an_i
         */
        RealMatrix3 _T_mat;

        
        MAST::FieldFunction<RealMatrixX>* _T_mat_function;
    };
    
}

#endif // __mast__local_elem_base__

