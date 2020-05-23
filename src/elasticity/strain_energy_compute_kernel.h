
#ifndef _mast_strain_energy_compute_kernel_h_
#define _mast_strain_energy_compute_kernel_h_

// MAST includes
#include "mesh/fe_compute_kernel.h"

namespace MAST {

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
    virtual inline value_type value() const override { }
    
protected:
    
    
    const MAST::ComputeKernel<material_value_type> *_material;
    const fe_var_type                              *_fe_var_data;
};

}

#endif // _mast_strain_energy_compute_kernel_h_
