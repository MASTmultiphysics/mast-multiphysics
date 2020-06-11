
#ifndef _mast_elasticity_elem_operations_h_
#define _mast_elasticity_elem_operations_h_


// MAST includes
#include "base/computation_base.h"
#include "mesh/quadrature.h"
#include "mesh/fe.h"
#include "mesh/fe_geometry_basis_derived_data.h"
#include "mesh/fe_var_data.h"
#include "elasticity/strain_energy_compute_kernel.h"

namespace MAST {

template <typename ScalarType, typename ContextType>
class ConstantSectionMaterial {
public:
    
    ConstantSectionMaterial(){}
    virtual ~ConstantSectionMaterial(){}
    virtual inline void value(const ContextType& c, typename EigenMatrix<ScalarType>::type& m) const {
        libmesh_assert_equal_to(m.rows(), 3);
        libmesh_assert_equal_to(m.cols(), 3);
        m(0,0) = 1.; m(0, 1) = .5;
        m(1,0) = .5; m(1, 1) = 1.;
        m(2,2) = 0.75;
    }
};


template <typename BasisScalarType, typename NodalScalarType, typename SolScalarType, typename ContextType>
struct ElasticityTraits {
    
    using context_type          = ContextType;
    using basis_scalar_type     = BasisScalarType;
    using nodal_scalar_type     = NodalScalarType;
    using sol_scalar_type       = SolScalarType;
    using var_scalar_type       = typename MAST::DeducedScalarType<NodalScalarType, SolScalarType>::type;
    using view_traits           = EigenTraits;
    using section_property_type = MAST::ConstantSectionMaterial<var_scalar_type, ContextType>;
    using quadrature_type       = MAST::libMeshQuadrature<ContextType>;
    using fe_basis_type         = MAST::libMeshFE<BasisScalarType, ContextType>;
    using fe_derived_type       = MAST::FEGeometryBasisDerivedData<BasisScalarType, NodalScalarType, view_traits, ContextType>;
    using fe_var_type           = MAST::FEVarData<basis_scalar_type, nodal_scalar_type, sol_scalar_type, view_traits, ContextType>;
};


struct ElasticityElemOperationsContext {
    
    ElasticityElemOperationsContext():
    qp      (0),
    elem    (nullptr),
    sys     (nullptr),
    physics (nullptr),
    sol     (nullptr) {}
    
    uint_type qp;
    const libMesh::Elem* elem;
    const MAST::NonlinearSystem* sys;
    const MAST::PhysicsDisciplineBase* physics;
    const libMesh::NumericVector<Real>* sol;

    inline Real nodal_coord(uint_type nd, uint_type d) const {
        libmesh_assert_msg(elem, "Element not initialized.");
        return elem->point(nd)(d);
    }
    
    inline uint_type n_nodes() const {
        libmesh_assert_msg(elem, "Element not initialized.");
        return elem->n_nodes();
    }
};


template <typename Traits>
class ElasticityElemOperations:
public MAST::ComputationBase<Traits, typename Traits::context_type> {
    
public:
    
    using context_type           = typename Traits::context_type;
    using var_scalar_type        = typename Traits::var_scalar_type;
    using quadrature_type        = typename Traits::quadrature_type;
    using fe_basis_type          = typename Traits::fe_basis_type;
    using fe_derived_type        = typename Traits::fe_derived_type;
    using fe_var_type            = typename Traits::fe_var_type;
    using section_property_type  = typename Traits::section_property_type;
    
    ElasticityElemOperations():
    MAST::ComputationBase<Traits, context_type>(),
    _initialized   (false),
    _sys           (nullptr),
    _physics       (nullptr),
    _X             (nullptr),
    _quadrature    (nullptr),
    _fe_basis      (nullptr),
    _fe_derived    (nullptr),
    _fe_var_data   (nullptr),
    _strain_energy (nullptr)
    { }
    
    virtual ~ElasticityElemOperations() { this->clear();}
    
    virtual void init(const MAST::NonlinearSystem& sys,
                      const MAST::PhysicsDisciplineBase& physics,
                      const libMesh::NumericVector<Real>& sol) {
        
        libmesh_assert(!_initialized);
        
        
        // setup the volume compute kernels
        _sys     = &sys;
        _physics = &physics;
        _X       = &sol;
        
        _quadrature    = new quadrature_type("quadrature");
        _fe_basis      = new fe_basis_type("fe_basis");
        _fe_basis->set_quadrature(*_quadrature);
        _fe_derived    = new fe_derived_type("fe_derived");
        _fe_derived->set_fe_basis(*_fe_basis);
        _fe_var_data   = new fe_var_type("u_vars");
        _fe_var_data->set_fe_shape_data(*_fe_derived);
        _strain_energy = new MAST::Element2D::LinearContinuumStrainEnergy<Traits, context_type>;
        _strain_energy->set_fe_var_data(*_fe_var_data);
        //_strain_energy->set_material(*_material);
        
        _quadrature->init_quadrature(libMesh::QGAUSS, 2);
        libMesh::FEType fe_type(1, libMesh::LAGRANGE);
        _fe_basis->init(_quadrature->get_libmesh_object(), fe_type);
        
        _initialized = true;
    }
    
    virtual void clear() {
        
        if (_initialized) {
         
            _sys     = nullptr;
            _physics = nullptr;
            delete    _quadrature;    _quadrature = nullptr;
            delete      _fe_basis;      _fe_basis = nullptr;
            delete    _fe_derived;    _fe_derived = nullptr;
            delete   _fe_var_data;   _fe_var_data = nullptr;
            delete _strain_energy; _strain_energy = nullptr;
        }
    }
    
    virtual void reinit(const libMesh::Elem& e) {
        
        libmesh_assert_msg(_initialized, "Object not initialized.");
        
        context_type c;
        c.elem    = &e;
        c.sys     = _sys;
        c.physics = _physics;
        c.sol     = _X;
        
        _fe_basis->reinit(e);
        _fe_derived->reinit(c);
        _fe_var_data->init(c);
    }
    
    virtual void compute(EigenVector<var_scalar_type>& res,
                         EigenMatrix<var_scalar_type>* jac) {}
    virtual void compute_sensitivity(const MAST::FunctionBase& f,
                                     EigenVector<var_scalar_type>& res) {}
    virtual void compute_complex_step_jacobian() {}
    virtual void compute_auto_diff_jacobian() {}
    
protected:
    
    bool  _initialized;
    const MAST::NonlinearSystem* _sys;
    const MAST::PhysicsDisciplineBase* _physics;
    const libMesh::NumericVector<Real>* _X;
    
    quadrature_type                      *_quadrature;
    fe_basis_type                        *_fe_basis;
    fe_derived_type                      *_fe_derived;
    fe_var_type                          *_fe_var_data;
    typename MAST::Element2D::LinearContinuumStrainEnergy<Traits, context_type>  *_strain_energy;
};
}

#endif // _mast_elasticity_elem_operations_h_
