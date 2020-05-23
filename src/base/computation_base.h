
#ifndef _mast_computation_base_h_
#define _mast_computation_base_h_

// MAST includes
#include "base/compute_kernel_base.h"

namespace MAST {

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

}

#endif // _mast_computation_base_h_
