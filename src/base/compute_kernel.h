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

class ComputeKernelDataBase {

public:

    ComputeKernelDataBase(const std::string& nm,
                          const std::vector<uint_type>& d):
    _name                 (nm),
    _dim                  (d)
    { }

    ComputeKernelDataBase(MAST::ComputeKernelDataBase& d):
    _name                 (d._name),
    _dim                  (d._dim)
    { }

    virtual ~ComputeKernelDataBase() {}
    
    inline uint_type n_dim() const { return _dim.size();}

protected:

    const std::string            _name;
    const std::vector<uint_type> _dim;
};



template <typename ScalarType>
class IndexableComputeKernelData: public MAST::ComputeKernelDataBase {
    
public:
    
    IndexableComputeKernelData(const std::string& nm,
                               const std::vector<uint_type>& d):
    MAST::ComputeKernelDataBase(nm, d),
    _view                (nullptr)
    { }
    
    IndexableComputeKernelData(MAST::IndexableComputeKernelData<ScalarType>& d):
    MAST::ComputeKernelDataBase(d) { }
    
    virtual ~IndexableComputeKernelData()
    { if (_view) delete _view; }


    inline void set_outer_dimensions(std::vector<uint_type> idx,
                                     uint_type buffer_size=0);
    
    template <typename ViewType> inline ViewType
    get_local_view(const std::vector<uint_type>& idx) {
        
        // make sure this object has been initialized
        libmesh_assert(_view);
        
        return _view->get_view_slice(idx);
    }
    
protected:

    MAST::View<ScalarType>       *_view;
    std::vector<uint_type>       _outer_idx_size;
};



/*! provides access to the element solution vector  through a memory view. */
template <typename SolScalarType>
class FieldVariable {
  
public:

    FieldVariable();
    virtual ~FieldVariable();
    /*! partial derivative of u with respect to time at */
    inline SolScalarType u() = 0;
    /*! partial derivative of u with respect to time at */
    inline SolScalarType du_dx(uint_type x_i) = 0;
    /*! partial derivative of u with respect to time at */
    inline SolScalarType du_dt() = 0;

};



// Forward decleration
template <typename> class ComputeKernelDerivative;


template <typename ValueType>
class ComputeKernel: public MAST::ComputeKernelBase {

public:
    
    using derivative_kernel_type = MAST::ComputeKernelDerivative<ValueType>;
    
    ComputeKernel(const std::string& nm):
    MAST::ComputeKernelBase(nm),
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
    
    virtual inline ValueType value() const = 0;

protected:
    
    derivative_kernel_type                  *_derivative_kernel;
};



template <typename ValueType>
class ComputeKernelDerivative: public MAST::ComputeKernelBase {

public:
    
    using primal_kernel_type   = MAST::ComputeKernel<ValueType>;

    ComputeKernelDerivative (const std::string& nm):
    MAST::ComputeKernelBase (nm),
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
    virtual inline ValueType value() const = 0;
    
protected:
    
    primal_kernel_type *_primal_kernel;
    const MAST::FunctionBase       *_f;
};

}

#endif /* __mast_compute_kernel_h__ */

