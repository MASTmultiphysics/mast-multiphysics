
#ifndef _mast_febasis_h_
#define _mast_febasis_h_

// MAST includes
#include "base/compute_kernel.h"
#include "mesh/quadrature.h"

namespace MAST {


template <typename ScalarType, typename Traits, typename ContextType>
class FEBasis: ComputeKernelBase<ContextType> {
    
public:

    using basis_scalar_type = ScalarType;
    
    FEBasis(const std::string& nm, const bool executable):
    MAST::ComputeKernelBase<ContextType>(nm, executable),
    _quadrature   (nullptr)
    { }
    virtual ~FEBasis() {}

    virtual inline void set_quadrature(const MAST::Quadrature<ScalarType, ContextType>& q) {
        
        libmesh_assert_msg(!_quadrature, "Quadrature already initialized.");
        _quadrature = &q;
    }
    
    virtual inline uint_type dim() const = 0;
    virtual inline uint_type order() const = 0;
    virtual inline uint_type n_basis() const = 0;
    
    virtual inline uint_type n_q_points() const {
        
        libmesh_assert_msg(_quadrature, "Quadrature not initialized.");
        return _quadrature->n_points();
    }
    
    virtual inline ScalarType qp_coord(uint_type qp, uint_type x_i) const
    {
        libmesh_assert_msg(_quadrature, "Quadrature not initialized.");
        return _quadrature->qp_coord(qp, x_i);
    }

    virtual inline ScalarType qp_weight(uint_type qp) const
    {
        libmesh_assert_msg(_quadrature, "Quadrature not initialized.");
        return _quadrature->weight(qp);
    }

    virtual inline ScalarType phi(uint_type qp, uint_type phi_i) const = 0;
    
    virtual inline ScalarType dphi_dxi(uint_type qp, uint_type phi_i, uint_type xi_i) const = 0;

protected:
    
    const MAST::Quadrature<ScalarType, ContextType> *_quadrature;
};


/*! serves as a wrapper around libMesh FEBase object to provide */
template <typename ScalarType, typename ContextType>
class libMeshFE: public FEBasis<ScalarType, EigenTraits, ContextType> {
    
public:
    
    using scalar_type       = ScalarType;
    using basis_scalar_type = typename MAST::FEBasis<ScalarType, EigenTraits, ContextType>::basis_scalar_type;
    
    libMeshFE(const std::string& nm):
    MAST::FEBasis<scalar_type, EigenTraits, ContextType>(nm, false),
    _compute_dphi_dxi   (false),
    _own_pointer        (false),
    _fe                 (nullptr)
    { }

    virtual ~libMeshFE() {
        
        if (_own_pointer)
            delete _fe;
    }
    
    inline void set_compute_dphi_dxi(bool f) { _compute_dphi_dxi = f;}
    
    virtual inline void init(libMesh::FEBase& fe) {

        libmesh_assert_msg(!_fe, "Object already initialized");
        
        _fe          = &fe;
        _own_pointer = false;

        const libMesh::FEMap& m = fe.get_fe_map();
        if (_compute_dphi_dxi)
            this->_dphi_dxi = {&m.get_dphidxi_map(), &m.get_dphideta_map(), &m.get_dphidzeta_map()};
        else
            _dphi_dxi.clear();
    }


    virtual inline void init(libMesh::QBase& q, const libMesh::FEType fe_type) {
        
        libmesh_assert_msg(!_fe, "Object already initialized");

        _fe = libMesh::FEBase::build(q.get_dim(), fe_type).release();
        _own_pointer = true;
        _fe->attach_quadrature_rule(&q);
    }

    
    virtual inline void clear() {

        if (_own_pointer)
            delete _fe;
        
        _own_pointer = false;
        _fe          = nullptr;
        _dphi_dxi.clear();
    }

    
    virtual inline void reinit(const libMesh::Elem& e) {
        
        libmesh_assert_msg(_fe, "Object not initialized");
        
        _fe->get_phi();
        if (_compute_dphi_dxi)
            _fe->get_JxW();
        
        _fe->reinit(&e);
        
        const libMesh::FEMap& m = _fe->get_fe_map();
        
        if (_compute_dphi_dxi)
            this->_dphi_dxi = {&m.get_dphidxi_map(), &m.get_dphideta_map(), &m.get_dphidzeta_map()};
        else
            _dphi_dxi.clear();
    }
    
    
    virtual inline uint_type dim() const override {
        
        libmesh_assert_msg(_fe, "FE not initialized.");

        return _fe->get_dim();
    }
    
    virtual inline uint_type order() const override {

        libmesh_assert_msg(_fe, "FE not initialized.");

        return _fe->get_order();
    }
    
    virtual inline uint_type n_basis() const override {

        libmesh_assert_msg(_fe, "FE not initialized.");
        
        return _fe->n_shape_functions();
    }
    
    virtual inline basis_scalar_type phi(uint_type qp, uint_type phi_i) const override
    {
        libmesh_assert_msg(_fe, "FE not initialized.");
        
        return _fe->get_phi()[phi_i][qp];
    }
    
    virtual inline basis_scalar_type
    dphi_dxi(uint_type qp, uint_type phi_i, uint_type xi_i) const override
    {
        libmesh_assert_msg(_compute_dphi_dxi, "FE not initialized with derivatives.");
        
        return (*_dphi_dxi[xi_i])[phi_i][qp];
    }
    
protected:
    
    bool _compute_dphi_dxi;
    bool _own_pointer;
    libMesh::FEBase* _fe;
    std::vector<const std::vector<std::vector<basis_scalar_type>>*> _dphi_dxi;
};




}

#endif // _mast_febasis_h_
