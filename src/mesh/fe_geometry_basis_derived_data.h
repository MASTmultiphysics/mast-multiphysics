
#ifndef _mast_fe_geometry_basis_derived_data_h_
#define _mast_fe_geometry_basis_derived_data_h_

// MAST includes
#include "mesh/fe_shape_data_base.h"


namespace MAST {

/*! This provides the derivative of shape functions when the FE basis also forms the basis for geometry interpolation used to interpolate
 nodal locations. Typically, Lagrange shape functoins are used for this purpose */
template <typename BasisScalarType, typename NodalScalarType, typename Traits, typename ContextType>
class FEGeometryBasisDerivedData:
public MAST::FEShapeDataBase<BasisScalarType, NodalScalarType, Traits, ContextType> {
    
public:
    
    FEGeometryBasisDerivedData(const std::string& nm):
    MAST::FEShapeDataBase<BasisScalarType, NodalScalarType, Traits, ContextType>(nm) {}
    virtual ~FEGeometryBasisDerivedData() {}
    
protected:

};



/*! This provides the derivative of shape functions when the FE basis also forms the basis for geometry interpolation used to interpolate
 nodal locations. Typically, Lagrange shape functoins are used for this purpose */
template <typename BasisScalarType, typename NodalScalarType, typename ContextType>
class FEGeometryBasisDerivedData<BasisScalarType, NodalScalarType, EigenTraits, ContextType>:
public MAST::FEShapeDataBase<BasisScalarType, NodalScalarType, EigenTraits, ContextType> {
    
public:
            
    FEGeometryBasisDerivedData(const std::string& nm):
    MAST::FEShapeDataBase<BasisScalarType, NodalScalarType, EigenTraits, ContextType>(nm) {}
    
    virtual ~FEGeometryBasisDerivedData() {}
    
    virtual inline void reinit(const ContextType& c) override {
        
        uint_type
        nq  = this->n_q_points();
        const int_type
        d   = this->spatial_dim();
        
        this->_xyz     = EigenMatrix<NodalScalarType>::type::Zero(nq, d);
        this->_detJ    = EigenVector<NodalScalarType>::type::Zero(nq);
        this->_detJxW  = EigenVector<NodalScalarType>::type::Zero(nq);
        this->_dx_dxi  = EigenMatrix<NodalScalarType>::type::Zero(nq, d*d);
        this->_dxi_dx  = EigenMatrix<NodalScalarType>::type::Zero(nq, d*d);
        this->_dphi_dx = EigenMatrix<NodalScalarType>::type::Zero(nq, d*d*this->n_basis());
    }

    virtual inline void reinit_for_side(const ContextType& c, uint_type s) override { }
    
    
protected:

};

}

#endif // _mast_fe_geometry_basis_derived_data_h_
