
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

template <typename ContextType>
class ComputeKernelBase {

public:

    ComputeKernelBase(const std::string& nm,
                      const bool executable):
    _nm          (nm),
    _executable  (executable) {}
    
    virtual ~ComputeKernelBase() {}
    
    virtual inline bool is_executable() const { return _executable;}
    virtual inline bool depends_on(const MAST::ComputeKernelBase<ContextType>& d) const
    { return _dependency.count(&d);}
    virtual inline const std::set<const MAST::ComputeKernelBase<ContextType>*>& get_dependencies() const
    { return _dependency;}

protected:

    virtual inline void _add_dependency(const MAST::ComputeKernelBase<ContextType>& d) { _dependency.insert(&d);}

    const std::string _nm;
    const bool        _executable;
    std::set<const MAST::ComputeKernelBase<ContextType>*> _dependency;
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
