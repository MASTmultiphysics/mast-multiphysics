
#ifndef _mast_fe_basis_derived_data_h_
#define _mast_fe_basis_derived_data_h_

// MAST includes
#include "mesh/fe_shape_data_base.h"

namespace MAST {

/*! This provides the derivative of shape functions when the FE basis is different from that used for interpolation of geometry. Typically,
 this is any non-Lagrange basis.*/
template <typename BasisScalarType, typename NodalScalarType, typename ViewTraits>
class FEBasisDerivedData: public MAST::FEShapeDataBase<BasisScalarType, NodalScalarType> {
  
public:
    
    using fe_geom_derived_type = MAST::FEGeometryBasisDerivedData<BasisScalarType, NodalScalarType, ViewTraits>;
    
    FEBasisDerivedData(const std::string& nm):
    MAST::FEShapeDataBase<BasisScalarType, NodalScalarType>(nm),
    _fe_geom(nullptr)
    {}
    
    virtual ~FEBasisDerivedData() {}
    
    virtual void inline set_fe_geom_derived_data(const fe_geom_derived_type& fe_geom)
    { _fe_geom = &fe_geom;}
    
    virtual inline NodalScalarType         xyz(uint_type qp, uint_type x_i) const override
    { return _fe_geom->xyz(qp, x_i);}
    
    virtual inline NodalScalarType        detJ(uint_type qp) const override
    { return _fe_geom->detJ(qp);}
    
    virtual inline NodalScalarType      detJxW(uint_type qp) const override
    { return _fe_geom->detJxW(qp);}

    virtual inline NodalScalarType      dx_dxi(uint_type qp, uint_type   x_i, uint_type xi_i) const override
    { return _fe_geom->dx_dxi(qp, x_i, xi_i);}

    virtual inline NodalScalarType      dxi_dx(uint_type qp, uint_type   x_i, uint_type xi_i) const override
    { return _fe_geom->dxi_dx(qp, x_i, xi_i);}

    virtual inline NodalScalarType     dphi_dx(uint_type qp, uint_type phi_i, uint_type x_i) const override
    { return _dphi_dx(qp, phi_i, x_i);}

protected:
    
    const fe_geom_derived_type             *_fe_geom;
    
    typename ViewTraits::dphi_dx_view_type _dphi_dx;
};



/*! This provides the derivative of shape functions when the FE basis is different from that used for interpolation of geometry. Typically,
 this is any non-Lagrange basis.*/
template <typename BasisScalarType, typename NodalScalarType>
class FEBasisDerivedData<BasisScalarType, NodalScalarType, MAST::EigenFEShapeDataViewTraits<NodalScalarType>>:
public MAST::FEShapeDataBase<BasisScalarType, NodalScalarType> {
    
public:
    
    using fe_derived_data_type = MAST::FEGeometryBasisDerivedData<BasisScalarType, NodalScalarType, MAST::EigenFEShapeDataViewTraits<NodalScalarType>>;
    
    FEBasisDerivedData(const std::string& nm):
    MAST::FEShapeDataBase<BasisScalarType, NodalScalarType>(nm),
    _fe_geom (nullptr)
    {}
    
    virtual ~FEBasisDerivedData() {}

    virtual void inline set_fe_geom_derived_data(const fe_derived_data_type& fe_geom)
    { _fe_geom = &fe_geom;}

    virtual inline NodalScalarType         xyz(uint_type qp, uint_type x_i) const override
    { return _fe_geom->xyz(qp, x_i);}
    
    virtual inline NodalScalarType        detJ(uint_type qp) const override
    { return _fe_geom->detJ(qp);}
    
    virtual inline NodalScalarType      detJxW(uint_type qp) const override
    { return _fe_geom->detJxW(qp);}

    virtual inline NodalScalarType      dx_dxi(uint_type qp, uint_type   x_i, uint_type xi_i) const override
    { return _fe_geom->dx_dxi(qp, x_i, xi_i);}

    virtual inline NodalScalarType      dxi_dx(uint_type qp, uint_type   x_i, uint_type xi_i) const override
    { return _fe_geom->dxi_dx(qp, x_i, xi_i);}

    virtual inline NodalScalarType     dphi_dx(uint_type qp, uint_type phi_i, uint_type x_i) const override
    { return _dphi_dx(qp, x_i*this->spatial_dim()+phi_i);}

protected:
    
    const fe_derived_data_type& _fe_geom;
    
    typename EigenFEShapeDataViewTraits<NodalScalarType>::dphi_dx_view_type _dphi_dx;
};
}

#endif // _mast_fe_basis_derived_data_h_
