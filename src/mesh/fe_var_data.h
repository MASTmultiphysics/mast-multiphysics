
#ifndef _mast_fe_var_data_h_
#define _mast_fe_var_data_h_

// MAST includes
#include "base/compute_kernel_base.h"


namespace MAST {

/*! provides access to the element solution vector  through a memory view. */
template <typename BasisScalarType, typename NodalScalarType, typename SolScalarType, typename Traits, typename ContextType>
class FEVarData: public MAST::ComputeKernelBase<ContextType> {
  
public:

    FEVarData(const std::string& nm): MAST::ComputeKernelBase<ContextType>(nm, false) {}
    virtual ~FEVarData() {}

};




template <typename BasisScalarType, typename NodalScalarType, typename SolScalarType, typename ContextType>
class FEVarData<BasisScalarType, NodalScalarType, SolScalarType, EigenTraits, ContextType>:
public MAST::ComputeKernelBase<ContextType> {
  
public:
    
    using var_scalar_type       = typename MAST::DeducedScalarType<NodalScalarType, SolScalarType>::type;
    using vector_type           = typename EigenTraits::vector_type<var_scalar_type>;
    using matrix_type           = typename EigenTraits::vector_type<var_scalar_type>;
    using fe_shape_data_type    = typename MAST::FEShapeDataBase<BasisScalarType, NodalScalarType, EigenTraits, ContextType>;
    
    FEVarData(const std::string& nm):
    MAST::ComputeKernelBase<ContextType>(nm, false),
    _compute_du_dx   (nullptr),
    _fe              (nullptr),
    _coeffs          (nullptr)
    {}
    virtual ~FEVarData() {}

    virtual inline void init() {
        
    }
    
    virtual inline void set_compute_du_dx(bool f) { _compute_du_dx = f;}
    virtual inline void set_fe_shape_data(const fe_shape_data_type& fe) { _fe = &fe;}
    virtual inline const fe_shape_data_type& get_fe_shape_data() const {
        
        libmesh_assert_msg(_fe, "FE pointer not initialized.");
        return *_fe;
    }
    virtual inline void set_fe_coefficient_view(const vector_type& coeffs) { _coeffs = &coeffs;}
    virtual inline uint_type n_q_points() const { return _u.size();}
    virtual inline var_scalar_type u(uint_type qp, uint_type comp) const  { return _u(qp, comp);}
    virtual inline var_scalar_type du_dx(uint_type qp, uint_type comp, uint_type x_i) const
    { return _du_dx(qp, _fe->spatial_dim()*comp + x_i);}

protected:
    
    bool                               _compute_du_dx;
    const fe_shape_data_type          *_fe;
    const vector_type                 *_coeffs;
    matrix_type                        _u;
    matrix_type                        _du_dx;
};

}

#endif // _mast_fe_var_data_h_
