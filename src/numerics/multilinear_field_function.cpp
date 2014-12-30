
// MAST includes
#include "numerics/multilinear_field_function.h"


template <typename ValType>
MAST::MultilinearInterpolation<ValType>::
MultilinearInterpolation(const std::string& nm,
                         MAST::FieldFunction<ValType>* abscissa,
                         std::map<Real, MAST::FieldFunction<ValType>*>& values):
MAST::FieldFunction<ValType>(nm),
_abscissa(abscissa),
_values(values) {
    
    // make sure that the size of the provided values is finite
    libmesh_assert(values.size() > 0);
    
    this->_functions.insert(_abscissa->master());
    
    typename std::map<Real, MAST::FieldFunction<ValType>*>::iterator
    it = values.begin(), end = values.end();
    
    // tell the function that it is dependent on the provided functions
    for ( ; it != end; it++)
        this->_functions.insert(it->second->master());
}



template <typename ValType>
MAST::MultilinearInterpolation<ValType>::
MultilinearInterpolation(const MAST::MultilinearInterpolation<ValType>& o):
MAST::FieldFunction<ValType>(o),
_abscissa(o._abscissa->clone().release()),
_values(o._values) {

    this->_functions.insert(_abscissa->master());
    
    typename std::map<Real, MAST::FieldFunction<ValType>*>::iterator
    it = _values.begin(), end = _values.end();
    
    // tell the function that it is dependent on the provided functions
    for ( ; it != end; it++)
        this->_functions.insert(it->second->master());
}



template <typename ValType>
MAST::MultilinearInterpolation<ValType>::~MultilinearInterpolation<ValType>() {
    
    delete _abscissa;
}



template <typename ValType>
std::auto_ptr<MAST::FieldFunction<ValType> >
MAST::MultilinearInterpolation<ValType>::clone() const {
    
    return std::auto_ptr<MAST::FieldFunction<ValType> >
    (new MAST::MultilinearInterpolation<ValType>(*this));
}


template <typename ValType>
void
MAST::MultilinearInterpolation<ValType>::
operator() (const libMesh::Point& p, Real t, Real& v) const {
    
    //
    // the following is used for calculation of the return value
    //   f(x) is defined for x for each x0 < x < x1
    //   if   x <= x0,      f(x) = f(x0)
    //   if   x0 < x < x1,  f(x) is interpolated
    //   if   x >= x1,      f(x) = f(x1)
    //
    
    typename std::map<Real, MAST::FieldFunction<ValType>*>::const_iterator
    it1, it2;
    typename std::map<Real, MAST::FieldFunction<ValType>*>::const_reverse_iterator
    rit = _values.rbegin();
    it1  = _values.begin();
    
    ValType x=0.;
    (*_abscissa)(p, t, x);
    
    // check the lower bound
    if (x <=  it1->first) {
        (*it1->second)(p, t, v);
    }
    // check the upper bound
    else if (x >=  rit->first) {
        (*rit->second)(p, t, v);
    }
    else {
        // if it gets here, the ordinate is in between the provided range
        it2 = _values.lower_bound(x);
        // this cannot be the first element of the map
        libmesh_assert(it2 != _values.begin());
        // it2 provides the upper bound. The lower bound is provided by the
        // preceding iterator
        it1 = it2--;
        Real f0 = 0., f1 = 0.;
        (*it1->second)(p, t, f0);
        (*it2->second)(p, t, f1);
        // now interpolate
        v =  (f0 +
              (x - it1->first)/(it2->first - it1->first) *
              (f1-f0));
    }
}




template <typename ValType>
void
MAST::MultilinearInterpolation<ValType>::
derivative(const MAST::DerivativeType d,
           const MAST::FunctionBase& f,
           const libMesh::Point& p,
           const Real t,
           ValType& v) const {
    
    //
    // the following is used for calculation of the return value
    //   f(x) is defined for x for each x0 < x < x1
    //   if   x <= x0,      f(x) = f(x0)
    //   if   x0 < x < x1,  f(x) is interpolated
    //   if   x >= x1,      f(x) = f(x1)
    //
    
    typename std::map<Real, MAST::FieldFunction<ValType>*>::const_iterator
    it1, it2;
    typename std::map<Real, MAST::FieldFunction<ValType>*>::const_reverse_iterator
    rit = _values.rbegin();
    it1  = _values.begin();
    
    ValType x=0., dx=0.;
    (*_abscissa)(p, t, x);

    // check the lower bound
    if (x <=  it1->first) {
        // since f(x) = f0 for x<x0, df/dx = 0
        v = 0.;
    }
    // check the upper bound
    else if (x >=  rit->first) {
        // since f(x) = fn for x>xn, df/dx = 0
        v = 0.;
    }
    else {
        // if it gets here, the ordinate is in between the provided range
        it2 = _values.lower_bound(x);
        // this cannot be the first element of the map
        libmesh_assert(it2 != _values.begin());
        // it2 provides the upper bound. The lower bound is provided by the
        // preceding iterator
        it1 = it2--;
        Real f0 = 0., f1 = 0.;
        it1->second->derivative(d, f, p, t, f0);
        it2->second->derivative(d, f, p, t, f1);
        // now interpolate: this calculates the term
        // dv/dp = sum_i df_i/dp eta
        // where eta is the interpolation
        v =  (f0 +
              (x - it1->first)/(it2->first - it1->first) *
              (f1-f0));
        
        // the other component of sensitivity is from dependence of x on p
        _abscissa->derivative(d, f, p, t, dx);
        (*it1->second)(p, t, f0);
        (*it2->second)(p, t, f1);
        v +=  (dx/(it2->first - it1->first) *
               (f1-f0));
    }
}



template class MAST::MultilinearInterpolation<Real>;
