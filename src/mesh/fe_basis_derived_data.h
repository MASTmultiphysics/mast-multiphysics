
#ifndef _mast_fe_basis_derived_data_h_
#define _mast_fe_basis_derived_data_h_

// MAST includes
#include "mesh/fe_shape_data_base.h"

namespace MAST {

/*! This provides the derivative of shape functions when the FE basis is different from that used for interpolation of geometry. Typically,
 this is any non-Lagrange basis.*/

template <typename BasisScalarType, typename NodalScalarType, typename Traits, typename ContextType>
class FEBasisDerivedData:
public MAST::FEShapeDataBase<BasisScalarType, NodalScalarType, Traits, ContextType> {
    
public:
    
    FEBasisDerivedData(const std::string& nm):
    MAST::FEShapeDataBase<BasisScalarType, NodalScalarType, Traits, ContextType>(nm, false),
    _fe_geom(nullptr)
    {}
    
    virtual ~FEBasisDerivedData() {}
    
};



/*! This provides the derivative of shape functions when the FE basis is different from that used for interpolation of geometry. Typically,
 this is any non-Lagrange basis.*/
template <typename BasisScalarType, typename NodalScalarType, typename ContextType>
class FEBasisDerivedData<BasisScalarType, NodalScalarType, EigenTraits>:
public MAST::FEShapeDataBase<BasisScalarType, NodalScalarType, EigenTraits, ContextType> {
    
public:
    
    using fe_derived_data_type = MAST::FEGeometryBasisDerivedData<BasisScalarType, NodalScalarType, EigenTraits>;
    
    FEBasisDerivedData(const std::string& nm):
    MAST::FEShapeDataBase<BasisScalarType, NodalScalarType, EigenTraits, ContextType>(nm),
    _fe_geom (nullptr)
    {}
    
    virtual ~FEBasisDerivedData() {}
    
    virtual void inline set_fe_geom_derived_data(const fe_derived_data_type& fe_geom)
    {
        libmesh_assert_msg(!_fe_geom, "FE object already initialized");
        
        _fe_geom = &fe_geom;
    }
    
    virtual inline NodalScalarType         xyz(uint_type qp, uint_type x_i) const override
    {
        libmesh_assert_msg(_fe_geom, "FE Geometry Basis not initialized.");
        
        return _fe_geom->xyz(qp, x_i);
    }
    
    virtual inline NodalScalarType        detJ(uint_type qp) const override
    {
        libmesh_assert_msg(_fe_geom, "FE Geometry Basis not initialized.");
        
        return _fe_geom->detJ(qp);
    }
    
    virtual inline NodalScalarType      detJxW(uint_type qp) const override
    {
        libmesh_assert_msg(_fe_geom, "FE Geometry Basis not initialized.");
        
        return _fe_geom->detJxW(qp);
    }
    
    virtual inline NodalScalarType      dx_dxi(uint_type qp, uint_type   x_i, uint_type xi_i) const override
    {
        libmesh_assert_msg(_fe_geom, "FE Geometry Basis not initialized.");
        
        return _fe_geom->dx_dxi(qp, x_i, xi_i);
    }
    
    virtual inline NodalScalarType      dxi_dx(uint_type qp, uint_type   x_i, uint_type xi_i) const override
    {
        libmesh_assert_msg(_fe_geom, "FE Geometry Basis not initialized.");
        
        return _fe_geom->dxi_dx(qp, x_i, xi_i);
    }
    
    
protected:
    
    const fe_derived_data_type& _fe_geom;
};
}

#endif // _mast_fe_basis_derived_data_h_
