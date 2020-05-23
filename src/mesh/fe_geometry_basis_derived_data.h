
#ifndef _mast_fe_geometry_basis_derived_data_h_
#define _mast_fe_geometry_basis_derived_data_h_

// MAST includes
#include "mesh/fe_shape_data_base.h"


namespace MAST {

/*! This provides the derivative of shape functions when the FE basis also forms the basis for geometry interpolation used to interpolate
 nodal locations. Typically, Lagrange shape functoins are used for this purpose */
template <typename BasisScalarType, typename NodalScalarType, typename ViewTraits>
class FEGeometryBasisDerivedData:
public MAST::FEShapeDataBase<BasisScalarType, NodalScalarType> {
    
public:
    
    FEGeometryBasisDerivedData(const std::string& nm):
    MAST::FEShapeDataBase<BasisScalarType, NodalScalarType>(nm) {}
    virtual ~FEGeometryBasisDerivedData() {}
    
    virtual inline NodalScalarType         xyz(uint_type qp, uint_type x_i) const override
    { return _xyz(qp, x_i);}
    
    virtual inline NodalScalarType        detJ(uint_type qp) const override
    { return _detJ(qp);}

    virtual inline NodalScalarType      detJxW(uint_type qp) const override
    { return _detJxW(qp);}

    virtual inline NodalScalarType      dx_dxi(uint_type qp, uint_type   x_i, uint_type xi_i) const override
    { return _dx_dxi(qp, x_i, xi_i);}

    virtual inline NodalScalarType      dxi_dx(uint_type qp, uint_type   x_i, uint_type xi_i) const override
    { return _dxi_dx(qp, x_i, xi_i);}

    virtual inline NodalScalarType     dphi_dx(uint_type qp, uint_type phi_i, uint_type x_i) const override
    { return _dphi_dx(qp, phi_i, x_i);}

protected:

    typename ViewTraits::xyz_view_type _xyz;
    typename ViewTraits::detJ_view_type _detJ;
    typename ViewTraits::detJxW_view_type _detJxW;
    typename ViewTraits::dx_dxi_view_type _dx_dxi;
    typename ViewTraits::dxi_dx_view_type _dxi_dx;
    typename ViewTraits::dphi_dx_view_type _dphi_dx;
};



/*! This provides the derivative of shape functions when the FE basis also forms the basis for geometry interpolation used to interpolate
 nodal locations. Typically, Lagrange shape functoins are used for this purpose */
template <typename BasisScalarType, typename NodalScalarType>
class FEGeometryBasisDerivedData<BasisScalarType, NodalScalarType, MAST::EigenFEShapeDataViewTraits<NodalScalarType>>:
public MAST::FEShapeDataBase<BasisScalarType, NodalScalarType> {
    
public:
    
    FEGeometryBasisDerivedData(const std::string& nm):
    MAST::FEShapeDataBase<BasisScalarType, NodalScalarType>(nm) {}
    
    virtual ~FEGeometryBasisDerivedData() {}
    
    virtual inline NodalScalarType         xyz(uint_type qp, uint_type x_i) const override
    { return _xyz(qp, x_i);}
    
    virtual inline NodalScalarType        detJ(uint_type qp) const override
    { return _detJ(qp);}

    virtual inline NodalScalarType      detJxW(uint_type qp) const override
    { return _detJxW(qp);}

    virtual inline NodalScalarType      dx_dxi(uint_type qp, uint_type   x_i, uint_type xi_i) const override
    { return _dx_dxi(qp, xi_i*this->spatial_dim()+x_i);}

    virtual inline NodalScalarType      dxi_dx(uint_type qp, uint_type   x_i, uint_type xi_i) const override
    { return _dxi_dx(qp, xi_i*this->spatial_dim()+x_i);}

    virtual inline NodalScalarType     dphi_dx(uint_type qp, uint_type phi_i, uint_type x_i) const override
    { return _dphi_dx(qp, x_i*this->spatial_dim()+phi_i);}

protected:

    typename EigenFEShapeDataViewTraits<NodalScalarType>::xyz_view_type _xyz;
    typename EigenFEShapeDataViewTraits<NodalScalarType>::detJ_view_type _detJ;
    typename EigenFEShapeDataViewTraits<NodalScalarType>::detJxW_view_type _detJxW;
    typename EigenFEShapeDataViewTraits<NodalScalarType>::dx_dxi_view_type _dx_dxi;
    typename EigenFEShapeDataViewTraits<NodalScalarType>::dxi_dx_view_type _dxi_dx;
    typename EigenFEShapeDataViewTraits<NodalScalarType>::dphi_dx_view_type _dphi_dx;
};

}

#endif // _mast_fe_geometry_basis_derived_data_h_
