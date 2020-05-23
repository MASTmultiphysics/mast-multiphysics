
#ifndef _mast_febasis_h_
#define _mast_febasis_h_

// MAST includes
#include "base/compute_kernel.h"
#include "mesh/quadrature.h"

namespace MAST {

template <typename ScalarType>
struct EigenFEShapeDataViewTraits {
    
    using xyz_view_type     = Eigen::Matrix<ScalarType, Dynamic, 3>;
    using detJ_view_type    = Eigen::Matrix<ScalarType, Dynamic, 1>;
    using detJxW_view_type  = Eigen::Matrix<ScalarType, Dynamic, 1>;
    using dx_dxi_view_type  = Eigen::Matrix<ScalarType, Dynamic, 9>;
    using dxi_dx_view_type  = Eigen::Matrix<ScalarType, Dynamic, 9>;
    using dphi_dx_view_type = Eigen::Matrix<ScalarType, Dynamic, Dynamic>;
};



template <typename ScalarType>
struct EigenFESolDataViewTraits {
    
    using coefficient_view_type = Eigen::Matrix<ScalarType, Dynamic, 1>;
    using u_view_type           = Eigen::Matrix<ScalarType, Dynamic, Dynamic>;
    using du_dx_view_type       = Eigen::Matrix<ScalarType, Dynamic, Dynamic>;
};




template <typename ScalarType>
class FEBasis: ComputeKernelBase {
    
public:

    using basis_scalar_type = ScalarType;
    
    FEBasis(const std::string& nm):
    MAST::ComputeKernelBase(nm),
    _quadrature   (nullptr)
    { }
    virtual ~FEBasis() {}

    virtual inline void set_quadrature(const MAST::Quadrature<ScalarType>& q) { _quadrature = &q;}
    virtual inline uint_type dim() const = 0;
    virtual inline uint_type order() const = 0;
    virtual inline uint_type n_basis() const = 0;
    virtual inline ScalarType qp_coord(uint_type qp, uint_type x_i) const
    { return _quadrature->qp_coord(qp, x_i); }
    virtual inline ScalarType phi(uint_type qp, uint_type phi_i) const = 0;
    virtual inline ScalarType dphi_dxi(uint_type qp, uint_type phi_i, uint_type xi_i) const = 0;
    
protected:
    
    const MAST::Quadrature<ScalarType> *_quadrature;
};


/*! serves as a wrapper around libMesh FEBase object to provide */
class libMeshFE: public FEBasis<Real> {
                                
public:
    
    using scalar_type = Real;
    
    libMeshFE(const std::string& nm):
    MAST::FEBasis<scalar_type>(nm),
    _compute_dphi_dxi   (false),
    _fe                 (nullptr)
    { }
    virtual ~libMeshFE() { }
    
    inline void set_compute_dphi_dxi(bool f) { _compute_dphi_dxi = f;}
    
    virtual inline void reinit(const libMesh::FEBase& fe) {
        _fe = &fe;
        const libMesh::FEMap& m = fe.get_fe_map();
        _dphi_dxi = {&m.get_dphidxi_map(), &m.get_dphideta_map(), &m.get_dphidzeta_map()};
    }

    virtual inline void execute() override { }
    
    virtual inline uint_type dim() const override { return _fe->get_dim();}
    
    virtual inline uint_type order() const override {return _fe->get_order();}
    
    virtual inline uint_type n_basis() const override { return _fe->n_shape_functions();}
    
    virtual inline basis_scalar_type phi(uint_type qp, uint_type phi_i) const override
    {return _fe->get_phi()[phi_i][qp];}
    
    virtual inline basis_scalar_type
    dphi_dxi(uint_type qp, uint_type phi_i, uint_type xi_i) const override
    { return (*_dphi_dxi[xi_i])[phi_i][qp];}

protected:
    
    bool _compute_dphi_dxi;
    const libMesh::FEBase* _fe;
    std::vector<const std::vector<std::vector<basis_scalar_type>>*> _dphi_dxi;
};




}

#endif // _mast_febasis_h_
