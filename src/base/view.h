
#ifndef _mast_view_h_
#define _mast_view_h_

// C++ includes
#include <map>

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





class ViewBase {
  
public:
    ViewBase(const std::string& nm, const std::vector<uint_type>& sizes):
    _name        (nm),
    _dim         (sizes.size()),
    _size        (sizes) {
        
        // this should be a vector or higher dimensional data
        libmesh_assert_greater(_dim, 0);
        
        _slice_size.resize(_dim+1, 1);
        
        for (uint_type i=_dim-1; i>=0; i++)
            _slice_size[i] = _slice_size[i+1]*_size[i];
    }
    
    virtual ~ViewBase() {}

    inline const std::string& name() const { return _name;}

    virtual inline void set_outer_sizes(const std::vector<uint_type>& d)
    {
        libmesh_assert_equal_to(_outside_size.size(), 0);
        _outside_size = d;
        _total_size.resize(d.size()+_size.size(), 0);
        for (uint_type i=0; i<d.size(); i++) _total_size[i] = d[i];
        for (uint_type i=0; i<_size.size(); i++) _total_size[i+d.size()] = _size[i];
    }

    inline uint_type array_index(const std::vector<uint_type>& indices) {

        libmesh_assert_equal_to(indices.size(), _dim);
        
        uint_type n = 0;
        for (uint_type i=0; i<_dim; i++) {
            libmesh_assert_less(indices[i], _size[i]);
            n += indices[i] * _slice_size[i+1];
        }

        libmesh_assert_less(n, _slice_size[0]);
        
        return n;
    }
    
    inline virtual void init() = 0;
    
protected:

    const std::string            _name;
    const uint_type              _dim;
    const std::vector<uint_type> _size;
    std::vector<uint_type>       _slice_size;
    std::vector<uint_type>       _outside_size;
    std::vector<uint_type>       _total_size;
};


template <typename ScalarType>
class View: public MAST::ViewBase {

public:
    
    View(const std::string& nm, const std::vector<uint_type>& sizes):
    MAST::ViewBase(nm, sizes),
    _array       (nullptr) {

    }

    virtual ~View() { delete _array;}

    virtual inline void init () override {

        libmesh_assert_msg(!_array, "View array already initialized.");
        
        _array = new ScalarType[_slice_size[0]];
        
        for (uint_type i=0; i<_slice_size[0]; i++)
            _array[i] = ScalarType{};
    }
    
    inline ScalarType& operator() (const std::vector<uint_type>& indices) {
        
        libmesh_assert_equal_to(indices.size(), _dim);
        return _array[array_index(indices)];
    }

    inline const ScalarType operator() (const std::vector<uint_type>& indices) const {
        
        libmesh_assert_equal_to(indices.size(), _dim);
        return _array[array_index(indices)];
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
    
    ScalarType                  *_array;
};


class ComputationData {
  
public:
    
    ComputationData(){}
    
    virtual ~ComputationData() {}
    
    inline void add_data(MAST::ViewBase& v) {

        const std::string& nm = v.name();
        
        libmesh_assert_msg(_data_map.find(nm) == _data_map.end(), "Duplicate data name: "+nm);
        
        _data_map[nm] = &v;
    }

     template <typename ScalarType> MAST::View<ScalarType>& get_data(const std::string& nm) {
        
        std::map<std::string, MAST::ViewBase*>::iterator
        it  = _data_map.find(nm),
        end = _data_map.end();
        
        libmesh_assert_msg( it != end, "Invalid data name: "+nm);
        return static_cast<MAST::View<ScalarType>&>(*it->second);
    }
    
    virtual inline void set_outer_sizes_and_initialize(const std::vector<uint_type>& d)
    {
        std::map<std::string, MAST::ViewBase*>::iterator
        it  = _data_map.begin(),
        end = _data_map.end();

        for ( ; it != end; it++) {
            it->second->set_outer_sizes(d);
            it->second->init();
        }
    }

protected:
    
    std::map<std::string, ViewBase*> _data_map;
};

}

#endif // _mast_view_h_
