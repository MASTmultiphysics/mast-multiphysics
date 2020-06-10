//
//  compute_kernel.h
//  compute_kernels
//
//  Created by Manav Bhatia on 5/1/20.
//  Copyright © 2020 Manav Bhatia. All rights reserved.
//

#ifndef _mast_compute_kernel_h_
#define _mast_compute_kernel_h_

// C++ includes
#include <vector>
#include <string>

// MAST includes
#include "base/compute_kernel_base.h"
#include "base/view.h"


namespace MAST {


template <typename ValueType>
struct IsScalarType {
    using type = std::false_type;
};


template <>
struct IsScalarType<Real> {
    using type = std::true_type;
};


template <>
struct IsScalarType<std::complex<Real>> {
    using type = std::true_type;
};



template <typename KernelValueType, typename IsIndexable, typename IsScalarType = typename MAST::IsScalarType<KernelValueType>::type>
struct KernelReturnType {

};


template <typename KernelValueType>
struct KernelReturnType<KernelValueType, std::false_type, std::false_type> {
    
    using type = KernelValueType;
};


template <typename KernelValueType>
struct KernelReturnType<KernelValueType, std::false_type, std::true_type> {
    
    using type = KernelValueType;
};



template <typename KernelValueType>
struct KernelReturnType<KernelValueType,  std::true_type, std::false_type> {
    
    using type = Eigen::Map<KernelValueType>;
};


template <typename KernelValueType>
struct KernelReturnType<KernelValueType,  std::true_type, std::true_type> {
    
    using type = KernelValueType&;
};


template <typename KernelType, typename KernelViewType = typename KernelType::view_type, typename ContextType>
KernelViewType build_kernel_view(KernelType& k, ContextType& c) { }

 


// Forward decleration
template <typename> class ComputeKernelDerivative;


template <typename ContextType>
class ComputeKernel: public MAST::ComputeKernelBase<ContextType> {

public:
    
    using derivative_kernel_type = MAST::ComputeKernelDerivative<ContextType>;
    
    ComputeKernel(const std::string& nm, const bool executable):
    MAST::ComputeKernelBase<ContextType>(nm, executable),
    _derivative_kernel   (nullptr)
    {}
    
    virtual ~ComputeKernel() {}
    
    virtual inline void set_derivative_kernel(derivative_kernel_type& d)
    {
        libmesh_assert_msg(!_derivative_kernel, "Derivative kernel already set.");
        _derivative_kernel = &d;
    }
    
    virtual inline const derivative_kernel_type& get_derivative_kernel() const
    {
        libmesh_assert_msg(_derivative_kernel, "Derivative kernel not set");
        return *_derivative_kernel;
    }
        
protected:
    
    derivative_kernel_type     *_derivative_kernel;
};



template <typename ContextType>
class ComputeKernelDerivative: public MAST::ComputeKernelBase<ContextType> {

public:
    
    using primal_kernel_type   = MAST::ComputeKernel<ContextType>;

    ComputeKernelDerivative (const std::string& nm, const bool executable):
    MAST::ComputeKernelBase<ContextType> (nm),
    _primal_kernel          (nullptr),
    _f                      (nullptr)
    {}
    
    virtual ~ComputeKernelDerivative() {}
    
    virtual inline void  set_primal_kernel(primal_kernel_type& k)
    {
        libmesh_assert_msg(!_primal_kernel, "Primal kernel already set.");
        _primal_kernel = &k;
    }

    virtual inline primal_kernel_type& get_primal_kernel()
    {
        libmesh_assert_msg(_primal_kernel, "Primal kernel not set.");
        return *_primal_kernel;
    }

    virtual inline void  set_derivative_paramter(const MAST::FunctionBase& f);
    
protected:
    
    primal_kernel_type *_primal_kernel;
    const MAST::FunctionBase       *_f;
};



}

#endif /* __mast_compute_kernel_h__ */

