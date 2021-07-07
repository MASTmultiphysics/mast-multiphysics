
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
    using matrix_type           = typename EigenTraits::matrix_type<var_scalar_type>;
    using fe_shape_data_type    = typename MAST::FEShapeDataBase<BasisScalarType, NodalScalarType, EigenTraits, ContextType>;
    
    FEVarData(const std::string& nm):
    MAST::ComputeKernelBase<ContextType>(nm, false),
    _compute_du_dx   (nullptr),
    _fe              (nullptr),
    _coeffs          (nullptr)
    {}
    virtual ~FEVarData() {}

    virtual inline void init(const ContextType& c,
                             const std::vector<uint_type>& active_vars,
                             const vector_type* coeffs = nullptr) {

        this->_init_coefficients(c, active_vars, coeffs);
        this->_init_variables(c);
    }

    virtual inline void init_for_complex_step_sens(const ContextType& c,
                                                   const uint_type complex_step_dof,
                                                   const std::vector<uint_type>& active_vars) {
        
        this->_init_coefficients(c, active_vars, nullptr);

        // the provided coefficients shoudl not already have
        libmesh_assert_equal_to(_coeffs->imag().norm(), 0.);
        
        // add complex-step perturbation to the complex-step dof specified
        libmesh_assert_less_equal(complex_step_dof, _coeffs->size());
        this->_add_complex_perturbation<var_scalar_type>(_coeff_vec, complex_step_dof);
        
        this->_init_variables(c);
    }

    inline void clear_coeffs_and_vars() {
        
        _coeffs = nullptr;
        _coeff_vec.setZero();
        
        _u.setZero();
        _du_dx.setZero();
    }
    
    virtual inline void set_compute_du_dx(bool f) { _compute_du_dx = f;}
    
    virtual inline void set_fe_shape_data(const fe_shape_data_type& fe) {
        
        libmesh_assert_msg(!_fe, "FE pointer already initialized.");
        _fe = &fe;
    }
    
    virtual inline const fe_shape_data_type& get_fe_shape_data() const {
        
        libmesh_assert_msg(_fe, "FE pointer not initialized.");
        return *_fe;
    }
    
    virtual inline uint_type n_q_points() const {
        
        libmesh_assert_msg(_fe, "FE pointer not initialized");
        return _fe->n_q_points();
    }
    
    virtual inline var_scalar_type u(uint_type qp, uint_type comp) const  {
        
        libmesh_assert_msg(_coeffs, "Object not initialized");
        libmesh_assert_less_equal(comp, _active_vars.size());

        return _u(qp, comp);
    }
    
    virtual inline var_scalar_type du_dx(uint_type qp, uint_type comp, uint_type x_i) const
    {
        libmesh_assert_msg(_compute_du_dx, "Object not initialized with du/dx");
        libmesh_assert_msg(_coeffs, "Object not initialized");
        libmesh_assert_less_equal(comp, _active_vars.size());
        
        return _du_dx(qp, _fe->spatial_dim()*comp + x_i);
    }

protected:
    
    inline void _init_coefficients(const ContextType& c,
                                   const std::vector<uint_type>& active_vars,
                                   const vector_type* coeffs = nullptr) {
        
        libmesh_assert_msg(!_coeffs, "Object already initialized");
        libmesh_assert_equal_to(_fe->spatial_dim(), c.elem->dim());
        libmesh_assert_less_equal(active_vars.size(), c.sys->n_vars());
        
        uint_type
        d        = _fe->spatial_dim(),
        n_qp     = _fe->n_q_points(),
        n_basis  = _fe->n_basis(),
        n_comp   = active_vars.size();

        _active_vars = active_vars;

        // first get the coefficients. If the vector coeffs is provided then
        // use that, else extract it from the context.
        if (coeffs) {
            
            _coeffs = coeffs;
            _coeff_vec.setZero();
        }
        else {
            
            std::vector<libMesh::dof_id_type> dofs;

            _coeff_vec = vector_type::Zero(n_comp*n_basis);

            for (uint_type i=0; i<n_comp; i++) {
                
                c.sys->get_dof_map().dof_indices(c.elem, dofs, active_vars[i]);
                
                for (uint_type j=0; j<dofs.size(); j++)
                    _coeff_vec(n_basis*i+j) = c.sol->el(dofs[j]);
            }
            
            _coeffs = &_coeff_vec;
        }

        libmesh_assert_equal_to(_coeffs->size(), n_comp*n_basis);
    }
    
    
    inline void _init_variables(const ContextType& c) {
        
        libmesh_assert_msg(_coeffs, "Coefficients not initialized");
        libmesh_assert_equal_to(_fe->spatial_dim(), c.elem->dim());
        libmesh_assert_less_equal(_active_vars.size(), c.sys->n_vars());
        
        uint_type
        d        = _fe->spatial_dim(),
        n_qp     = _fe->n_q_points(),
        n_basis  = _fe->n_basis(),
        n_comp   = _active_vars.size();

        _u = matrix_type::Zero(n_qp, n_comp);
        
        // now, initialize the solution value and derivatives.
        for (uint_type i=0; i<n_qp; i++)
            for (uint_type j=0; j<n_comp; j++)
                for (uint_type k=0; k<n_basis; k++)
                    _u(i, j) += _fe->phi(i, k) * (*_coeffs)(j*n_basis+k);
        
        if (_compute_du_dx) {
            
            _du_dx = matrix_type::Zero(n_qp, n_comp*d);

            for (uint_type i=0; i<n_qp; i++)
                for (uint_type j=0; j<n_comp; j++)
                    for (uint_type l=0; l<d; l++)
                        for (uint_type k=0; k<n_basis; k++)
                            _du_dx(i, d*j+l) += _fe->dphi_dx(i, k, l) * (*_coeffs)(j*n_basis+k);

        }
        else
            _du_dx.setZero();
    }
    
    template <typename ValType>
    inline void _add_complex_perturbation(vector_type& v, uint_type i) {
        
        libmesh_assert_less_equal(i, v.size());
        
        v(i) += Complex(0, ComplexStepDelta);
    }

    template <>
    inline void _add_complex_perturbation<Real>(vector_type& v, uint_type i) {

        libmesh_assert_msg(false, "Complex perturbation cannot be added to real vectors");
    }

    bool                               _compute_du_dx;
    const fe_shape_data_type          *_fe;
    const vector_type                 *_coeffs;
    std::vector<uint_type>             _active_vars;
    vector_type                        _coeff_vec;
    matrix_type                        _u;
    matrix_type                        _du_dx;
};

}

#endif // _mast_fe_var_data_h_
