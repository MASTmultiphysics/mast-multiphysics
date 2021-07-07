
#ifndef _mast_fe_shape_data_base_h_
#define _mast_fe_shape_data_base_h_

// MAST includes
#include "mesh/fe.h"

namespace MAST {

template <typename BasisScalarType, typename NodalScalarType, typename Traits, typename ContextType>
class FEShapeDataBase: public MAST::ComputeKernelBase<ContextType> {

public:
    
    FEShapeDataBase(const std::string& nm):
    MAST::ComputeKernelBase<ContextType>(nm, false) {}
    
    virtual ~FEShapeDataBase() {}

};


    
template <typename BasisScalarType, typename NodalScalarType, typename ContextType>
class FEShapeDataBase<BasisScalarType, NodalScalarType, EigenTraits, ContextType>:
public MAST::ComputeKernelBase<ContextType> {
    
public:

    using dphi_dx_type = typename Eigen::Map<const typename EigenVector<NodalScalarType>::type>;
    
    FEShapeDataBase(const std::string& nm):
    MAST::ComputeKernelBase<ContextType>(nm, false),
    _compute_xyz       (false),
    _compute_Jac       (false),
    _compute_Jac_inv   (false),
    _compute_detJ      (false),
    _compute_JxW       (false),
    _compute_dphi_dx   (false),
    _compute_normal    (false),
    _spatial_dim       (0),
    _fe_basis          (nullptr)
    { }
    
    virtual ~FEShapeDataBase() {}
    
    inline void      set_compute_xyz(bool f) { _compute_xyz = f;}
    
    inline void      set_compute_Jac(bool f) {
        
        _compute_Jac = f;
        if (f) this->set_compute_xyz(true);
    }

    inline void      set_compute_Jac_inverse(bool f) {
        
        _compute_Jac_inv = f;
        if (f) this->set_compute_Jac(true);
    }

    inline void     set_compute_detJ(bool f) {
        
        _compute_detJ = f;
        if (f) this->set_compute_Jac(true);
    }
    
    inline void   set_compute_detJxW(bool f) {
        
        _compute_JxW = f;
        if (f) this->set_compute_detJ(true);
    }
    
    inline void  set_compute_dphi_dx(bool f) {
        
        _compute_dphi_dx = f;
        if (f) this->set_compute_Jac_inverse(true);
    }
    
    inline void   set_compute_normal(bool f) {
        
        _compute_normal = f;
        if (f) this->set_compute_Jac(true);
    }
    
    inline void  set_fe_basis(MAST::FEBasis<BasisScalarType, EigenTraits, ContextType>& basis)
    {
        libmesh_assert_msg(!_fe_basis, "FE Basis already initialized.");
        
        _fe_basis = &basis;
    }
    
    virtual inline void reinit(const ContextType& c) = 0;
    virtual inline void reinit_for_side(const ContextType& c, uint_type s) = 0;

    virtual inline uint_type             ref_dim() const {
        
        libmesh_assert_msg(_fe_basis, "FE Basis not initialized.");
        return _fe_basis->dim();
    }
    
    virtual inline uint_type         spatial_dim() const {
        
        libmesh_assert_msg(_fe_basis, "FE Basis not initialized.");
        return _spatial_dim;
    }
    
    virtual inline uint_type               order() const {
    
        libmesh_assert_msg(_fe_basis, "FE Basis not initialized.");
        return _fe_basis->order();
    }
    
    virtual inline uint_type             n_basis() const {
    
        libmesh_assert_msg(_fe_basis, "FE Basis not initialized.");
        return _fe_basis->n_basis();
    }
    
    virtual inline BasisScalarType         phi(uint_type qp, uint_type phi_i) const
    {
        libmesh_assert_msg(_fe_basis, "FE Basis not initialized.");
         return _fe_basis->phi(qp, phi_i);
    }
    
    virtual inline BasisScalarType    dphi_dxi(uint_type qp, uint_type phi_i, uint_type xi_i) const
    {
    
        libmesh_assert_msg(_fe_basis, "FE Basis not initialized.");

        return _fe_basis->dphi_dxi(qp, phi_i, xi_i);
    }
    
    virtual inline uint_type        n_q_points() const { return _fe_basis->n_q_points();}

    virtual inline NodalScalarType  node_coord(uint_type nd, uint_type x_i) const
    {
        libmesh_assert_msg(_compute_xyz, "Nodal and QPoint locations not requested");
        return _node_coord(nd, x_i);
    }

    virtual inline NodalScalarType         xyz(uint_type qp, uint_type x_i) const
    {
        libmesh_assert_msg(_compute_xyz, "Nodal and QPoint locations not requested");
        return _xyz(qp, x_i);
    }
    
    virtual inline NodalScalarType        detJ(uint_type qp) const
    {
        libmesh_assert_msg(_compute_detJ, "Jacobian computation not requested");
        return _detJ(qp);
    }

    virtual inline NodalScalarType      detJxW(uint_type qp) const
    {
        libmesh_assert_msg(_compute_JxW, "JxW computation not requested");
        return _detJxW(qp);
    }

    virtual inline NodalScalarType      dx_dxi(uint_type qp, uint_type   x_i, uint_type xi_i) const
    {
        libmesh_assert_msg(_compute_Jac, "Jacobian computation not requested");
        return _dx_dxi(qp, xi_i*this->spatial_dim()+x_i);
    }

    virtual inline NodalScalarType      dxi_dx(uint_type qp, uint_type   x_i, uint_type xi_i) const
    {
        libmesh_assert_msg(_compute_Jac_inv, "Jacobian inverse computation not requested");
        return _dxi_dx(qp, xi_i*this->spatial_dim()+x_i);
    }

    inline const dphi_dx_type
    dphi_dx(uint_type qp, uint_type x_i) const
    {
        libmesh_assert_msg(_compute_dphi_dx, "Jacobian inverse computation not requested");

        return dphi_dx_type(_dphi_dx.row(qp).segment(x_i*this->n_basis(), this->n_basis()).data(),
                            this->n_basis());
    }
    
    virtual inline NodalScalarType     dphi_dx(uint_type qp, uint_type phi_i, uint_type x_i) const
    {
        libmesh_assert_msg(_compute_dphi_dx, "Jacobian inverse computation not requested");

        return _dphi_dx(qp, x_i*this->n_basis()+phi_i);
    }

protected:

    bool _compute_xyz;
    bool _compute_Jac;
    bool _compute_Jac_inv;
    bool _compute_detJ;
    bool _compute_JxW;
    bool _compute_dphi_dx;
    bool _compute_normal;
    uint_type _spatial_dim;
    
    MAST::FEBasis<BasisScalarType, EigenTraits, ContextType> *_fe_basis;
    
    typename EigenMatrix<NodalScalarType>::type       _node_coord;
    typename EigenMatrix<NodalScalarType>::type       _xyz;
    typename EigenVector<NodalScalarType>::type       _detJ;
    typename EigenVector<NodalScalarType>::type       _detJxW;
    typename EigenMatrix<NodalScalarType>::type       _dx_dxi;
    typename EigenMatrix<NodalScalarType>::type       _dxi_dx;
    typename EigenMatrix<NodalScalarType>::type       _dphi_dx;
};

}
#endif // _mast_fe_shape_data_base_h_
