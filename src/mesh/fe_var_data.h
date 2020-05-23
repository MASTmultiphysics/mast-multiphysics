
#ifndef _mast_fe_var_data_h_
#define _mast_fe_var_data_h_

// MAST includes
#include "base/compute_kernel_base.h"


namespace MAST {

/*! provides access to the element solution vector  through a memory view. */
template <typename BasisScalarType, typename NodalScalarType, typename SolScalarType, typename ViewTraits, typename ContextType>
class FEVarData: public MAST::ComputeKernelBase<ContextType> {
  
public:

    using var_scalar_type       = typename MAST::DeducedScalarType<NodalScalarType, SolScalarType>::type;
    using coefficient_view_type = typename ViewTraits::coefficient_view_type;
    using fe_shape_data_type    = typename MAST::FEShapeDataBase<BasisScalarType, NodalScalarType, ContextType>;

    FEVarData(const std::string& nm): MAST::ComputeKernelBase<ContextType>(nm) {}
    virtual ~FEVarData() {}
    virtual inline void execute(ContextType& c) override {}
    inline void set_compute_du_dx(bool f) { _compute_du_dx = f;}
    inline void
    set_fe_shape_data(const fe_shape_data_type& fe) { _fe = &fe;}
    inline void set_fe_coefficient_view(const coefficient_view_type& coeffs) { _coeffs = coeffs;}
    inline var_scalar_type u(uint_type qp, uint_type comp) { return _u(qp, comp);}
    inline var_scalar_type du_dx(uint_type qp, uint_type comp, uint_type x_i) { return _du_dx(qp, comp, x_i);}

protected:
    
    bool                                                *_compute_du_dx;
    const fe_shape_data_type                            *_fe;
    const typename ViewTraits::coefficient_view_type    _coeffs;
    typename ViewTraits::u_view_type                    _u;
    typename ViewTraits::du_dx_view_type                _du_dx;
};




template <typename BasisScalarType, typename NodalScalarType, typename SolScalarType, typename ContextType>
class FEVarData<BasisScalarType, SolScalarType, NodalScalarType, MAST::EigenFESolDataViewTraits<SolScalarType>, ContextType>:
public MAST::ComputeKernelBase<ContextType> {
  
public:
    
    using var_scalar_type       = typename MAST::DeducedScalarType<NodalScalarType, SolScalarType>::type;
    using coefficient_view_type = typename EigenFESolDataViewTraits<SolScalarType>::coefficient_view_type;
    using u_view_type           = typename MAST::EigenFESolDataViewTraits<var_scalar_type>::u_view_type;
    using du_dx_view_type       = typename MAST::EigenFESolDataViewTraits<var_scalar_type>::du_dx_view_type;
    using fe_shape_data_type    = typename MAST::FEShapeDataBase<BasisScalarType, NodalScalarType, ContextType>;
    
    FEVarData(const std::string& nm):
    MAST::ComputeKernelBase<ContextType>(nm),
    _compute_du_dx   (nullptr),
    _fe              (nullptr),
    _coeffs          (nullptr)
    {}
    virtual ~FEVarData() {}

    virtual inline void execute(ContextType& c) override {}
    inline void set_compute_du_dx(bool f) { _compute_du_dx = f;}
    inline void set_fe_shape_data(const fe_shape_data_type& fe) { _fe = &fe;}
    inline void set_fe_coefficient_view(const coefficient_view_type& coeffs) { _coeffs = &coeffs;}
    inline var_scalar_type u(uint_type qp, uint_type comp) { return _u(qp, comp);}
    inline var_scalar_type du_dx(uint_type qp, uint_type comp, uint_type x_i)
    { return _du_dx(qp, _fe->spatial_dim()*comp + x_i);}

protected:
    
    bool                               _compute_du_dx;
    const fe_shape_data_type          *_fe;
    const coefficient_view_type       *_coeffs;
    u_view_type                        _u;
    du_dx_view_type                    _du_dx;
};

}

#endif // _mast_fe_var_data_h_
