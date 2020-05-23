

#ifndef _mast_hookean_material_h_
#define _mast_hookean_material_h_

// MAST includes
#include "mesh/fe_compute_kernel.h"

namespace MAST {

template <typename Traits>
class HookeanMaterial:
public ComputeKernel<typename Traits::material_value_type, Traits> {

public:

    using value_type          = typename Traits::material_value_type;

    HookeanMaterial(const std::string& nm):
    MAST::ComputeKernel<typename Traits::material_value_type, Traits>(nm)
    { }
    
    virtual ~HookeanMaterial() { }
    
protected:
    
};

}

#endif // _mast_hookean_material_h_
