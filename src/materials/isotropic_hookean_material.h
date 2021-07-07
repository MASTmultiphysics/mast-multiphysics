

#ifndef _mast_isotropic_hookean_material_h_
#define _mast_isotropic_hookean_material_h_

// MAST includes
#include "mesh/fe_compute_kernel.h"

namespace MAST {

template <typename Traits, typename ContextType>
class IsotropicHookeanMaterial:
public ComputeKernel<typename Traits::material_value_type, Traits, ContextType> {

public:

    using scalar_type    = typename Traits::var_scalar_type;
    using value_type     = typename Traits::material_value_type;

    IsotropicHookeanMaterial(const std::string& nm):
    MAST::ComputeKernel<typename Traits::material_value_type, Traits>(nm),
    _E    (nullptr),
    _nu   (nullptr)
    { }
    
    virtual ~IsotropicHookeanMaterial() { }
    
    inline void set_properties(const MAST::ComputeKernel<scalar_type>& E,
                               const MAST::ComputeKernel<scalar_type>& nu) {
        
        libmesh_assert_msg(!_E && !_nu, "Property kernels should be cleared before setting.");
        
        _E   = &E;
        _nu  = &nu;
    }


    void inline operator(const ContextType& c, value_type& m) {
        
        libmesh_assert_msg(_E && _nu, "Property kernels should be provided before execute().");

        const scalar_type
        E  = _E->value(c),
        nu = _nu->value(c),
        G  = E/2./(1.+nu);

        m(0, 0) = m(1, 1) = m(2, 2) = E*(1.-nu)/(1.-nu-2.*nu*nu);
        m(0, 1) = m(0, 2) = m(1, 0) = m(1, 2) = m(2, 0) = m(2, 1) = E*nu/(1.-nu-2.*nu*nu);
        m(3, 3) = m(4, 4) = m(5, 5) = E/2./(1.+nu);
    }

    inline void execute_derivative(const ContextType& c, value_type& m) {
        
        libmesh_assert_msg(_E && _nu, "Property kernels should be provided before execute().");

        const scalar_type
        E  = _E->value(c),
        nu = _nu->value(c),
        G  = E/2./(1.+nu);

        MAST::init_view(_matrix, 6, 6, true);

        m(0, 0) = m(1, 1) = m(2, 2) = E*(1.-nu)/(1.-nu-2.*nu*nu);
        m(0, 1) = m(0, 2) = m(1, 0) = m(1, 2) = m(2, 0) = m(2, 1) = E*nu/(1.-nu-2.*nu*nu);
        m(3, 3) = m(4, 4) = m(5, 5) = E/2./(1.+nu);
    }

protected:

    const MAST::ComputeKernel<scalar_type> *_E;
    const MAST::ComputeKernel<scalar_type> *_nu;
};

}

#endif // _mast_isotropic_hookean_material_h_

