
#ifndef _mast_view_h_
#define _mast_view_h_

// MAST includes
#include "base/mast_data_types.h"



namespace MAST {


template <typename ScalarType, uint_type Dim>
struct EigenViewType {

    using eigen_matrix_type = Eigen::Matrix<ScalarType, Dynamic, Dynamic>;
    using view_type         = Eigen::Map<eigen_matrix_type>;
    using view_m1_type      = typename MAST::EigenViewType<ScalarType, Dim-1>::view_type;
    
    static inline view_type
    create_view(uint_type* dim, ScalarType* array) {
        uint_type
        n1 = dim[0],
        n2 = dim[1];
        for (unsigned int i=2; i<Dim; i++) n2 *= dim[i];
        
        return view_type(array, n1, n2);
    }

    /*! create slice of view at row idx */
    template <uint_type Level> static inline
    typename MAST::EigenViewType<ScalarType, Dim-Level>::view_type
    create_view_slice(const std::vector<uint_type> idx, uint_type* dim, ScalarType* array);

    template <> static inline
    typename MAST::EigenViewType<ScalarType, Dim-1>::view_type
    create_view_slice<1>(const std::vector<uint_type> idx, uint_type* dim, ScalarType* array)
    {
        
        libmesh_assert_equal_to(idx.size(), 1);
        libmesh_assert_less(idx[0], dim[0]);
        
        uint_type
        n2 = dim[1];
        
        for (unsigned int i=2; i<Dim; i++) n2 *= dim[i];

        // create a Dim-1 dimensional array starting at an offset
        return MAST::EigenViewType<ScalarType, Dim-1>::create_view(dim+1, array+idx[0]*n2);
    }
};


template <typename ScalarType>
struct EigenViewType<ScalarType, 1> {

    using eigen_matrix_type = Eigen::Matrix<ScalarType, Dynamic, 1>;
    using view_type         = Eigen::Map<eigen_matrix_type>;
    using view_m1_type      = ScalarType;

    template <uint_type Level> static inline ScalarType
    create_view_slice(const std::vector<uint_type> idx, uint_type* dim, ScalarType* array);

    template <> static inline ScalarType
    create_view_slice<1>(const std::vector<uint_type> idx, uint_type* dim, ScalarType* array) {

        libmesh_assert_less(idx.size(), 1);
        libmesh_assert_less(idx[0], dim[0]);
        
        // Dim-1 dimensional data is a scalar
        return array[idx[0]];
    }

};




template <typename ScalarType>
class View {

public:
    
    View(const std::string& nm, uint_type dim, uint_type n1, ...):
    _name        (nm),
    _dim         (dim),
    _dim_size    (nullptr),
    _array       (nullptr) {
        
        _dim_size   = new uint_type[dim];

        va_list args;
        va_start (args, n1);

        uint_type
        n = 1;
        
        for (uint_type i=0; i<dim; i++) {
            
            _dim_size[i] = va_arg(args, uint_type);
            n *= _dim[i];
        }

        _array = new ScalarType[n];
        
        va_end(args);
            
        for (uint_type i=0; i<n; i++)
            _array[i] = ScalarType{};
    }

    virtual ~View()
    {
        delete _dim_size;
        delete _array;
    }
    
    
    template <typename ViewType> inline ViewType get_view()
    { return MAST::EigenViewType<ScalarType, _dim>::create_view(_dim, _array);}

    template <typename ViewType> inline const ViewType get_view() const
    { return MAST::EigenViewType<ScalarType, _dim>::create_view(_dim, _array);}

    template <typename ViewType> inline ViewType get_view_slice(const std::vector<uint_type>& idx)
    { return MAST::EigenViewType<ScalarType, _dim>::create_view_slice<idx.size()>(idx, _dim, _array);}

    template <typename ViewType> inline const ViewType get_view_slice(const std::vector<uint_type>& idx) const
    { return MAST::EigenViewType<ScalarType, _dim>::create_view_slice<idx.size()>(idx, _dim, _array);}

    
protected:
    
    std::string       _name;
    uint_type         _dim;
    uint_type        *_dim_size;
    ScalarType       *_array;
};



}

#endif // _mast_view_h_
