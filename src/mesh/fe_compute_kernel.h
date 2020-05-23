
#ifndef _mast_fe_compute_kernel_h_
#define _mast_fe_compute_kernel_h_

// MAST includes
#include "base/compute_kernel.h"

namespace MAST {

template <typename ValueType, typename Traits, typename ContextType>
class FEComputeKernel: public MAST::ComputeKernel<ValueType, ContextType> {

public:

    using basis_scalar_type = typename Traits::basis_scalar_type;
    using nodal_scalar_type = typename Traits::nodal_scalar_type;
    using sol_scalar_type   = typename Traits::sol_scalar_type;
    
    FEComputeKernel(const std::string& nm):
    MAST::ComputeKernel<ValueType, ContextType>(nm) { }
    virtual ~FEComputeKernel() { }
};

}

#endif // _mast_fe_compute_kernel_h_
