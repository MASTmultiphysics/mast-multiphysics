

#ifndef _mast_quadrature_h_
#define _mast_quadrature_h_

// MAST includes
#include "base/compute_kernel.h"


namespace MAST {

template <typename ScalarType>
class Quadrature: public MAST::ComputeKernelBase {
    
public:

    using scalar_type = ScalarType;
    
    Quadrature(const std::string& nm): MAST::ComputeKernelBase(nm) {}
    virtual ~Quadrature() {}
    virtual inline uint_type dim() const = 0;
    virtual inline uint_type order() const = 0;
    virtual inline uint_type n_points() const = 0;
    virtual inline ScalarType qp_coord(uint_type qp, uint_type xi_i) const = 0;
    virtual inline ScalarType weight(uint_type qp) const = 0;
    
};



/*! This provides a  */
template <typename ScalarType, typename ViewTraits>
class MappedQuadrature: public MAST::Quadrature<ScalarType> {
    
public:

    using scalar_type = ScalarType;
    
    MappedQuadrature(const std::string& nm): MAST::Quadrature<ScalarType>(nm), _dim(0) {}
    virtual ~MappedQuadrature() {}
    virtual void set_data(const uint_type dim,
                          const typename ViewTraits::quadrature_point_type points,
                          const typename ViewTraits::quadrature_point_type weights) {
        _points = points;
        _weights = weights;
    }
    virtual inline uint_type dim() const override { return _dim;}
    virtual inline uint_type order() const override { /* should not be called*/ libmesh_assert(false);}
    virtual inline uint_type n_points() const override { return _weights.size();}
    virtual inline ScalarType qp_coord(uint_type qp, uint_type xi_i) const override { return _points(qp, xi_i);}
    virtual inline ScalarType weight(uint_type qp) const { return _weights(qp);}
    
protected:

    uint_type _dim;
    typename ViewTraits::quadrature_point_type _points;
    typename ViewTraits::quadrature_point_type _weights;
};



/*! serves as a wrapper around libMesh */
class libMeshQuadrature: public MAST::Quadrature<Real> {
    
public:
    
    /*! the quadrature object \p q  is expected to be initialized outside of this class. */
    libMeshQuadrature(const std::string& nm):
    MAST::Quadrature<Real>(nm),
    _q  (nullptr)
    {}
    virtual ~libMeshQuadrature() {}
    virtual inline void execute() override { }
    const libMesh::QBase& get_libmesh_object() const { return *_q;}
    virtual inline uint_type dim() const override { return _q->get_dim();}
    virtual inline uint_type order() const override { return _q->get_order();}
    virtual inline uint_type n_points() const override { return _q->n_points();}
    virtual inline scalar_type qp_coord(uint_type qp, uint_type xi_i) const override
    { return _q->get_points()[qp](xi_i);}
    virtual inline scalar_type weight(uint_type qp) const override { return _q->w(qp);}
    
protected:
    
    libMesh::QBase* _q;
};
}

#endif  //_mast_quadrature_h_
