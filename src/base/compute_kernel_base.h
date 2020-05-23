
#ifndef _mast_compute_kernel_base_h_
#define _mast_compute_kernel_base_h_

// C++ includes
#include <set>
#include <memory>

// MAST includes
#include "base/mast_data_types.h"


namespace MAST {

// Forward declerations
class FunctionBase;

class ComputeKernelBase {

public:

    ComputeKernelBase(const std::string& nm): _nm(nm) {}
    virtual ~ComputeKernelBase() {}
    virtual inline bool depends_on(const MAST::ComputeKernelBase& d) const { return _dependency.count(&d);}
    virtual inline const std::set<const MAST::ComputeKernelBase*>& get_dependencies() const { return _dependency;}
    virtual inline void init() {}
    virtual inline void pre_execute() {}
    virtual inline void post_execute() {}
    virtual inline void execute() = 0;

protected:

    virtual inline void _add_dependency(const MAST::ComputeKernelBase& d) { _dependency.insert(&d);}

    const std::string _nm;
    std::set<const MAST::ComputeKernelBase*> _dependency;
};


struct CurrentComputation {
    
    CurrentComputation(): time (0.), dt(0.), param(nullptr) {}
    
    std::vector<uint_type>  current_indices;
    Real                    time;
    Real                    dt;
    /*! parameter for which sensitivity is being computed */
    FunctionBase     *param;
};

}

#endif // _mast_compute_kernel_base_h_
