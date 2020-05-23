
#ifndef _mast_fe_shape_data_base_h_
#define _mast_fe_shape_data_base_h_

// MAST includes
#include "mesh/fe.h"

namespace MAST {

template <typename BasisScalarType, typename NodalScalarType>
class FEShapeDataBase: public MAST::ComputeKernelBase {
    
public:
    
    FEShapeDataBase(const std::string& nm):
    MAST::ComputeKernelBase(nm),
    _compute_xyz       (false),
    _compute_Jac       (false),
    _compute_detJ      (false),
    _compute_JxW       (false),
    _compute_dphi_dx   (false),
    _compute_normal    (false),
    _spatial_dim       (0),
    _fe_basis          (nullptr)
    {}
    virtual ~FEShapeDataBase() {}
    inline void      set_compute_xyz(bool f) { _compute_xyz = f;}
    inline void      set_compute_Jac(bool f) { _compute_Jac = f;}
    inline void     set_compute_detJ(bool f) { _compute_detJ = f;}
    inline void   set_compute_detJxW(bool f) { _compute_JxW = f;}
    inline void  set_compute_dphi_dx(bool f) { _compute_dphi_dx = f;}
    inline void   set_compute_normal(bool f) { _compute_normal = f;}
    inline void  set_fe_basis(MAST::FEBasis<BasisScalarType>& basis) { _fe_basis = &basis;}
    virtual inline void execute() override { }
    //virtual inline void reinit_for_side(MAST::FEBasis<BasisScalarType>& basis, uint_type s) = 0;
    virtual inline uint_type             ref_dim() const { return _fe_basis->dim();}
    virtual inline uint_type         spatial_dim() const { return _spatial_dim;}
    virtual inline uint_type               order() const { return _fe_basis->order();}
    virtual inline uint_type             n_basis() const { return _fe_basis->n_basis();}
    virtual inline BasisScalarType         phi(uint_type qp, uint_type phi_i) const
    { return _fe_basis->phi(qp, phi_i);}
    virtual inline BasisScalarType    dphi_dxi(uint_type qp, uint_type phi_i, uint_type xi_i) const
    { return _fe_basis->dphi_dxi(qp, phi_i, xi_i);}
    virtual inline NodalScalarType         xyz(uint_type qp, uint_type x_i) const = 0;
    virtual inline NodalScalarType        detJ(uint_type qp) const = 0;
    virtual inline NodalScalarType      detJxW(uint_type qp) const = 0;
    virtual inline NodalScalarType      dx_dxi(uint_type qp, uint_type   x_i, uint_type xi_i) const = 0;
    virtual inline NodalScalarType      dxi_dx(uint_type qp, uint_type   x_i, uint_type xi_i) const = 0;
    virtual inline NodalScalarType     dphi_dx(uint_type qp, uint_type phi_i, uint_type xi_i) const = 0;

protected:

    bool _compute_xyz;
    bool _compute_Jac;
    bool _compute_detJ;
    bool _compute_JxW;
    bool _compute_dphi_dx;
    bool _compute_normal;
    uint_type _spatial_dim;
    
    MAST::FEBasis<BasisScalarType> *_fe_basis;
};

}
#endif // _mast_fe_shape_data_base_h_
