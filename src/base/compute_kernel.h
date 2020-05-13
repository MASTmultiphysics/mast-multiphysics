//
//  compute_kernel.h
//  compute_kernels
//
//  Created by Manav Bhatia on 5/1/20.
//  Copyright © 2020 Manav Bhatia. All rights reserved.
//

#ifndef _mast_compute_kernel_h_
#define _mast_compute_kernel_h_

// C++ includes
#include <set>
#include <memory>

// MAST includes
#include "base/mast_data_types.h"

// libMesh includes
#include "libmesh/quadrature.h"
#include "libmesh/fe_base.h"
#include "libmesh/elem.h"

typedef unsigned int uint_type;

namespace MAST {


// Forward declerations
class FunctionBase;


class ComputeKernelBase {

public:

    ComputeKernelBase(const std::string& nm): _nm(nm) {}
    virtual ~ComputeKernelBase() {}
    virtual inline bool depends_on(const MAST::ComputeKernelBase& d) const { return _dependency.count(&d);}
    virtual inline const std::set<const MAST::ComputeKernelBase*>& get_dependencies() const { return _dependency;}
    virtual inline void pre_execute() {}
    virtual inline void post_execute() {}
    virtual inline void execute() = 0;

protected:

    virtual inline void _add_dependency(const MAST::ComputeKernelBase& d) { _dependency.insert(&d);}

    const std::string _nm;
    std::set<const MAST::ComputeKernelBase*> _dependency;
};


struct CurrentComputation {
    
    CurrentComputation(): qp(0), time (0.), dt(0.), param(nullptr) {}
    
    uint_type        qp;
    Real             time;
    Real             dt;
    /*! parameter for which sensitivity is being computed */
    FunctionBase     *param;
};


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


template <typename NodalScalarType, typename SolScalarType>
struct DeducedScalarType { };

template <>
struct DeducedScalarType<Real, Real> { using type = Real;};

template <>
struct DeducedScalarType<std::complex<Real>, Real> { using type = std::complex<Real>;};

template <>
struct DeducedScalarType<Real, std::complex<Real>> { using type = std::complex<Real>;};

template <>
struct DeducedScalarType<std::complex<Real>, std::complex<Real>> { using type = std::complex<Real>;};



template <typename ScalarType>
struct EigenFEShapeDataViewTraits {
    
    using xyz_view_type     = Eigen::Matrix<ScalarType, Dynamic, 3>;
    using detJ_view_type    = Eigen::Matrix<ScalarType, Dynamic, 1>;
    using detJxW_view_type  = Eigen::Matrix<ScalarType, Dynamic, 1>;
    using dx_dxi_view_type  = Eigen::Matrix<ScalarType, Dynamic, 9>;
    using dxi_dx_view_type  = Eigen::Matrix<ScalarType, Dynamic, 9>;
    using dphi_dx_view_type = Eigen::Matrix<ScalarType, Dynamic, Dynamic>;
    
    template <typename ViewType>
    void init_view(ViewType& view, uint_type n1);

    template <typename ViewType>
    void init_view(ViewType& view, uint_type n1, uint_type n2);
};



template <typename ScalarType>
struct EigenFESolDataViewTraits {
    
    using coefficient_view_type = Eigen::Matrix<ScalarType, Dynamic, 1>;
    using u_view_type           = Eigen::Matrix<ScalarType, Dynamic, Dynamic>;
    using du_dx_view_type       = Eigen::Matrix<ScalarType, Dynamic, Dynamic>;

    template <typename ViewType>
    void init_view(ViewType& view, uint_type n1);

    template <typename ViewType>
    void init_view(ViewType& view, uint_type n1, uint_type n2);
};



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



/*! provides access to the element solution vector  through a memory view. */
template <typename SolScalarType>
class FieldVariable {
  
public:

    FieldVariable();
    virtual ~FieldVariable();
    /*! partial derivative of u with respect to time at */
    inline SolScalarType u() = 0;
    /*! partial derivative of u with respect to time at */
    inline SolScalarType du_dx(uint_type x_i) = 0;
    /*! partial derivative of u with respect to time at */
    inline SolScalarType du_dt() = 0;

protected:
    
};



/*! provides access to the element solution vector  through a memory view. */
template <typename BasisScalarType, typename NodalScalarType, typename SolScalarType, typename ViewTraits>
class FEVarData: public MAST::ComputeKernelBase {
  
public:

    using var_scalar_type       = typename MAST::DeducedScalarType<NodalScalarType, SolScalarType>::type;
    using coefficient_view_type = typename ViewTraits::coefficient_view_type;

    FEVarData(const std::string& nm): MAST::ComputeKernelBase(nm) {}
    virtual ~FEVarData() {}
    virtual inline void execute() override {}
    inline void set_compute_du_dx(bool f) { _compute_du_dx = f;}
    inline void set_fe_shape_data(const MAST::FEShapeDataBase<BasisScalarType, NodalScalarType>& fe) { _fe = &fe;}
    inline void set_fe_coefficient_view(const coefficient_view_type& coeffs) { _coeffs = coeffs;}
    inline var_scalar_type u(uint_type qp, uint_type comp) { return _u(qp, comp);}
    inline var_scalar_type du_dx(uint_type qp, uint_type comp, uint_type x_i) { return _du_dx(qp, comp, x_i);}

protected:
    
    bool                                                             *_compute_du_dx;
    const MAST::FEShapeDataBase<BasisScalarType, NodalScalarType>    *_fe;
    const typename ViewTraits::coefficient_view_type                 _coeffs;
    typename ViewTraits::u_view_type                                 _u;
    typename ViewTraits::du_dx_view_type                             _du_dx;
};




template <typename BasisScalarType, typename NodalScalarType, typename SolScalarType>
class FEVarData<BasisScalarType, SolScalarType, NodalScalarType, MAST::EigenFESolDataViewTraits<SolScalarType>>:
public MAST::ComputeKernelBase {
  
public:
    
    using var_scalar_type       = typename MAST::DeducedScalarType<NodalScalarType, SolScalarType>::type;
    using coefficient_view_type = typename EigenFESolDataViewTraits<SolScalarType>::coefficient_view_type;
    using u_view_type           = typename MAST::EigenFESolDataViewTraits<var_scalar_type>::u_view_type;
    using du_dx_view_type       = typename MAST::EigenFESolDataViewTraits<var_scalar_type>::du_dx_view_type;

    
    FEVarData(const std::string& nm):
    MAST::ComputeKernelBase(nm),
    _compute_du_dx   (nullptr),
    _fe              (nullptr),
    _coeffs          (nullptr)
    {}
    virtual ~FEVarData() {}

    virtual inline void execute() override {}
    inline void set_compute_du_dx(bool f) { _compute_du_dx = f;}
    inline void set_fe_shape_data(const MAST::FEShapeDataBase<BasisScalarType, NodalScalarType>& fe) { _fe = &fe;}
    inline void set_fe_coefficient_view(const coefficient_view_type& coeffs) { _coeffs = &coeffs;}
    inline var_scalar_type u(uint_type qp, uint_type comp) { return _u(qp, comp);}
    inline var_scalar_type du_dx(uint_type qp, uint_type comp, uint_type x_i)
    { return _du_dx(qp, _fe->spatial_dim()*comp + x_i);}

protected:
    
    bool                                                             _compute_du_dx;
    const MAST::FEShapeDataBase<BasisScalarType, NodalScalarType>    *_fe;
    const coefficient_view_type                                      *_coeffs;
    u_view_type                                                      _u;
    du_dx_view_type                                                  _du_dx;
};



template <typename ValueType>
class ComputeKernel: public MAST::ComputeKernelBase {

public:
    
    ComputeKernel(const std::string& nm): MAST::ComputeKernelBase(nm) {}
    virtual ~ComputeKernel() {}
    virtual inline const ValueType& value() const = 0;
};



template <typename ValueType>
class ComputeKernelDerivative: public MAST::ComputeKernelBase {

public:
    
    ComputeKernelDerivative(const std::string& nm):
    MAST::ComputeKernelBase(nm),
    _f  (nullptr)
    {}
    virtual ~ComputeKernelDerivative() {}
    virtual inline void  set_derivative_paramter(const MAST::FunctionBase& f);
    virtual inline const ValueType& value() const = 0;
    
protected:
    
    const MAST::FunctionBase* _f;
};



template <typename ValueType, typename Traits>
class FEComputeKernel: public MAST::ComputeKernel<ValueType> {

public:

    using basis_scalar_type = typename Traits::basis_scalar_type;
    using nodal_scalar_type = typename Traits::nodal_scalar_type;
    using sol_scalar_type   = typename Traits::sol_scalar_type;
    
    FEComputeKernel(const std::string& nm): MAST::ComputeKernel<ValueType>(nm) { }
    virtual ~FEComputeKernel() { }
    virtual inline const ValueType& value() const = 0;
    //virtual inline const ValueType& derivative(const MAST::FunctionBase& f) const = 0;
};




/*! Collection of compute kernels with specified dependencies and data views*/
template <typename Traits>
class ComputationBase: public MAST::CurrentComputation {
    
public:
    
    using basis_scalar_type = typename Traits::basis_scalar_type;
    using nodal_scalar_type = typename Traits::nodal_scalar_type;
    using sol_scalar_type   = typename Traits::sol_scalar_type;

    ComputationBase(): MAST::CurrentComputation() { }
    virtual ~ComputationBase() { }
    void add_compute_kernel(MAST::ComputeKernelBase& c);
    /*! parses through all the compute kernels and their views to prepare a graph of dependency. */
    void prepare();
    void print_graph();
    void execute();
    
protected:
    
};


template <typename BasisScalarType, typename NodalScalarType, typename SolScalarType>
struct ElasticityComputeKernelTraits {

    using basis_scalar_type = BasisScalarType;
    using nodal_scalar_type = NodalScalarType;
    using sol_scalar_type   = SolScalarType;
    using var_scalar_type   = typename MAST::DeducedScalarType<NodalScalarType, SolScalarType>;
    using res_vector_type   = Eigen::Matrix<var_scalar_type, Dynamic, 1>;
    using jac_matrix_type   = Eigen::Matrix<var_scalar_type, Dynamic, Dynamic>;
    using material_val_type = Eigen::Matrix<var_scalar_type, Dynamic, Dynamic>;
    using quadrature_type   = MAST::libMeshQuadrature;
    using fe_basis_type     = MAST::libMeshFE;
    using fe_derived_traits = MAST::EigenFEShapeDataViewTraits<NodalScalarType>;
    using fe_derived_type   = MAST::FEGeometryBasisDerivedData<BasisScalarType, NodalScalarType, fe_derived_traits>;
    using fe_var_traits     = MAST::EigenFESolDataViewTraits<SolScalarType>;
    using fe_var_type       = typename MAST::FEVarData<basis_scalar_type, nodal_scalar_type, sol_scalar_type, fe_var_traits>;
};



template <typename Traits>
class StrainEnergy:
public FEComputeKernel<typename Traits::res_vector_type, Traits> {

public:

    using value_type          = typename Traits::res_vector_type;
    using basis_scalar_type   = typename Traits::basis_scalar_type;
    using nodal_scalar_type   = typename Traits::nodal_scalar_type;
    using sol_scalar_type     = typename Traits::sol_scalar_type;
    using fe_var_traits       = typename Traits::fe_var_traits;
    using fe_var_type         = typename Traits::fe_var_type;
    using material_value_type = typename Traits::material_val_type;

    StrainEnergy():
    MAST::FEComputeKernel<value_type, Traits>("strain_energy"),
    _material    (nullptr),
    _fe_var_data (nullptr)
    { }
    
    virtual ~StrainEnergy() { }
    virtual inline void set_material(const MAST::ComputeKernel<material_value_type>& material)
    { _material = &material;}
    virtual inline void set_fe_var_data(const fe_var_type& fe_data)
    { _fe_var_data = &fe_data;}
    virtual inline void execute() override { }
    virtual inline const value_type& value() const override { }
    
protected:
    
    
    const MAST::ComputeKernel<material_value_type> *_material;
    const fe_var_type                              *_fe_var_data;
};




template <typename Traits>
class ElasticityElemOperations: public MAST::ComputationBase<Traits> {

public:

    using value_type          = typename Traits::res_vector_type;
    using basis_scalar_type   = typename Traits::basis_scalar_type;
    using nodal_scalar_type   = typename Traits::nodal_scalar_type;
    using sol_scalar_type     = typename Traits::sol_scalar_type;
    using quadrature_type     = typename Traits::quadrature_type;
    using fe_view_traits      = typename Traits::fe_var_traits;
    using fe_basis_type       = typename Traits::fe_basis_type;
    using fe_derived_type     = typename Traits::fe_derived_type;
    using fe_var_type         = typename Traits::fe_var_type;
    using material_value_type = typename Traits::material_val_type;

    ElasticityElemOperations():
    MAST::ComputationBase<Traits>(),
    _initialized   (false),
    _sys           (nullptr),
    _physics       (nullptr),
    _X             (nullptr),
    _fe_var_data   (nullptr),
    _strain_energy (nullptr)
    { }
    virtual ~ElasticityElemOperations() {}
    virtual void init(const MAST::NonlinearSystem& sys,
                      const MAST::PhysicsDisciplineBase& physics,
                      const libMesh::NumericVector<Real>& X);
    virtual void clear() {}
    virtual void reinit(const libMesh::Elem& e) {}
    virtual void compute_residual() {}
    virtual void compute_jacobian() {}
    virtual void compute_complex_step_jacobian() {}
    virtual void compute_auto_diff_jacobian() {}
    
protected:
    
    bool  _initialized;
    const MAST::NonlinearSystem* _sys;
    const MAST::PhysicsDisciplineBase* _physics;
    const libMesh::NumericVector<Real>* _X;
    
    quadrature_type                      *_quadrature;
    fe_basis_type                        *_fe_basis;
    fe_derived_type                      *_fe_derived;
    fe_var_type                          *_fe_var_data;
    typename MAST::StrainEnergy<Traits>  *_strain_energy;
};


template <typename Traits>
void
MAST::ElasticityElemOperations<Traits>::init(const MAST::NonlinearSystem& sys,
                                             const MAST::PhysicsDisciplineBase& physics,
                                             const libMesh::NumericVector<Real>& X) {
    
    libmesh_assert(!_initialized);

    
    // setup the volume compute kernels
    _sys     = &sys;
    _physics = &physics;
    _X       = &X;
    
    
    _quadrature    = new quadrature_type("quadrature");
    _fe_basis      = new fe_basis_type("fe_basis");
    _fe_basis->set_quadrature(*_quadrature);
    _fe_derived    = new fe_derived_type("fe_derived");
    _fe_derived->set_fe_basis(*_fe_basis);
    _fe_var_data   = new fe_var_type("u_vars");
    _fe_var_data->set_fe_shape_data(*_fe_derived);
    _strain_energy = new MAST::StrainEnergy<Traits>;
    _strain_energy->set_fe_var_data(*_fe_var_data);
    //_strain_energy->set_material(*_material);
}


}
#endif /* __mast_compute_kernel_h__ */

