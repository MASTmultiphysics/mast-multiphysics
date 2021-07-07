
#ifndef _mast_strain_energy_compute_kernel_h_
#define _mast_strain_energy_compute_kernel_h_

// MAST includes
#include "mesh/fe_compute_kernel.h"
#include "numerics/fem_operator_matrix.h"


namespace MAST {

namespace Element2D {

template <typename NodalScalarType, typename VarScalarType, typename FEVarType>
void linear_continuum_strain(const FEVarType& fe_var,
                             const uint_type qp,
                             typename EigenVector<VarScalarType>::type& epsilon,
                             MAST::FEMOperatorMatrixBase<NodalScalarType>& Bmat_lin) {
    
    epsilon.setZero();
    
    const typename FEVarType::fe_shape_data_type
    &fe = fe_var.get_fe_shape_data();
    
    // make sure all matrices are the right size
    libmesh_assert_equal_to(epsilon.size(), 3);
    libmesh_assert_equal_to(Bmat_lin.m(), 3);
    libmesh_assert_equal_to(Bmat_lin.n(), 2*fe.n_basis());
    
    
    // linear strain operator
    Bmat_lin.set_shape_function(0, 0, fe.dphi_dx(qp, 0)); //  epsilon_xx = du/dx
    Bmat_lin.set_shape_function(2, 1, fe.dphi_dx(qp, 0)); //  gamma_xy = dv/dx + ...
        
    // linear strain operator
    Bmat_lin.set_shape_function(1, 1, fe.dphi_dx(qp, 1)); //  epsilon_yy = dv/dy
    Bmat_lin.set_shape_function(2, 0, fe.dphi_dx(qp, 1)); //  gamma_xy = du/dy + ...
    
    epsilon(0) = fe_var.du_dx(qp, 0, 0);  // du/dx
    epsilon(1) = fe_var.du_dx(qp, 1, 1);  // dv/dy
    epsilon(2) = fe_var.du_dx(qp, 0, 1) + fe_var.du_dx(qp, 1, 0);  // du/dy + dv/dx
}


template <typename ScalarType, typename FEVarType, typename ContextType>
void green_lagrange_strain_operator(const FEVarType& fe_var,
                                    const ContextType& qp,
                                    bool if_nonlinear,
                                    EigenMatrix<ScalarType>& E,
                                    EigenMatrix<ScalarType>& F,
                                    EigenVector<ScalarType>& epsilon,
                                    EigenMatrix<ScalarType>& mat_x,
                                    EigenMatrix<ScalarType>& mat_y,
                                    MAST::FEMOperatorMatrixBase<ScalarType>& Bmat_lin,
                                    MAST::FEMOperatorMatrixBase<ScalarType>& Bmat_nl_x,
                                    MAST::FEMOperatorMatrixBase<ScalarType>& Bmat_nl_y,
                                    MAST::FEMOperatorMatrixBase<ScalarType>& Bmat_nl_u,
                                    MAST::FEMOperatorMatrixBase<ScalarType>& Bmat_nl_v) {
    
    
    epsilon.setZero();
    mat_x.setZero();
    mat_y.setZero();

    const typename FEVarType::fe_shape_data_type
    &fe = fe_var.get_fe_shape_object();

    // make sure all matrices are the right size
    libmesh_assert_equal_to(epsilon.size(), 3);
    libmesh_assert_equal_to(mat_x.rows(), 3);
    libmesh_assert_equal_to(mat_x.cols(), 2);
    libmesh_assert_equal_to(mat_y.rows(), 3);
    libmesh_assert_equal_to(mat_y.cols(), 2);
    libmesh_assert_equal_to(Bmat_lin.m(), 3);
    libmesh_assert_equal_to(Bmat_lin.n(), 2*fe.n_basis());
    libmesh_assert_equal_to(Bmat_nl_x.m(), 2);
    libmesh_assert_equal_to(Bmat_nl_x.n(), 2*fe.n_basis());
    libmesh_assert_equal_to(Bmat_nl_y.m(), 2);
    libmesh_assert_equal_to(Bmat_nl_y.n(), 2*fe.n_basis());
    libmesh_assert_equal_to(F.cols(), 2);
    libmesh_assert_equal_to(F.rows(), 2);
    libmesh_assert_equal_to(E.cols(), 2);
    libmesh_assert_equal_to(E.rows(), 2);
    
    
    // now set the shape function values
    Bmat_lin.set_shape_function(0, 0, fe.dphi_dx(qp, 0)); //  epsilon_xx = du/dx
    Bmat_lin.set_shape_function(2, 1, fe.dphi_dx(qp, 0)); //  gamma_xy = dv/dx + ...
    
    // nonlinear strain operator in x
    Bmat_nl_x.set_shape_function(0, 0, fe.dphi_dx(qp, 0)); // du/dx
    Bmat_nl_x.set_shape_function(1, 1, fe.dphi_dx(qp, 0)); // dv/dx
    
    // nonlinear strain operator in u
    Bmat_nl_u.set_shape_function(0, 0, fe.dphi_dx(qp, 0)); // du/dx
    Bmat_nl_v.set_shape_function(0, 1, fe.dphi_dx(qp, 0)); // dv/dx
    
    // dN/dy
    Bmat_lin.set_shape_function(1, 1, fe.dphi_dx(qp, 1)); //  epsilon_yy = dv/dy
    Bmat_lin.set_shape_function(2, 0, fe.dphi_dx(qp, 1)); //  gamma_xy = du/dy + ...
    
    // nonlinear strain operator in y
    Bmat_nl_y.set_shape_function(0, 0, fe.dphi_dx(qp, 1)); // du/dy
    Bmat_nl_y.set_shape_function(1, 1, fe.dphi_dx(qp, 1)); // dv/dy
    
    // nonlinear strain operator in v
    Bmat_nl_u.set_shape_function(1, 0, fe.dphi_dx(qp, 1)); // du/dy
    Bmat_nl_v.set_shape_function(1, 1, fe.dphi_dx(qp, 1)); // dv/dy
    
    // prepare the deformation gradient matrix
    F.row(0) = fe_var.du_dx(qp, 0);
    F.row(1) = fe_var.du_dx(qp, 1);
    
    // this calculates the Green-Lagrange strain in the reference config
    E = 0.5*(F + F.transpose() + F.transpose() * F);
    
    // now, add this to the strain vector
    epsilon(0) = E(0,0);
    epsilon(1) = E(1,1);
    epsilon(2) = E(0,1) + E(1,0);
    
    // now initialize the matrices with strain components
    // that multiply the Bmat_nl terms
    mat_x(0, 0) =     fe_var.du_dx(qp, 0, 0);
    mat_x(0, 1) =     fe_var.du_dx(qp, 1, 0);
    mat_x(2, 0) =     fe_var.du_dx(qp, 0, 1);
    mat_x(2, 1) =     fe_var.du_dx(qp, 1, 1);
    
    mat_y(1, 0) =     fe_var.du_dx(qp, 0, 1);
    mat_y(1, 1) =     fe_var.du_dx(qp, 1, 1);
    mat_y(2, 0) =     fe_var.du_dx(qp, 0, 0);
    mat_y(2, 1) =     fe_var.du_dx(qp, 1, 0);
}


template <typename Traits, typename ContextType>
class LinearContinuumStrainEnergy:
public ComputeKernel<ContextType> {
    
public:
    
    using nodal_scalar_type     = typename Traits::nodal_scalar_type;
    using var_scalar_type       = typename Traits::var_scalar_type;
    using vector_type           = typename Traits::view_traits::template vector_type<var_scalar_type>;
    using matrix_type           = typename Traits::view_traits::template matrix_type<var_scalar_type>;
    using fe_var_type           = typename Traits::fe_var_type;
    using section_property_type = typename Traits::section_property_type;
    
    LinearContinuumStrainEnergy():
    MAST::ComputeKernel<ContextType>("strain_energy", false),
    _property    (nullptr),
    _fe_var_data (nullptr)
    { }
    
    virtual ~LinearContinuumStrainEnergy() { }

    virtual inline void
    set_section_property(const section_property_type& p) {
        
        libmesh_assert_msg(!_property, "Property already initialized.");
        
        _property = &p;
    }

    virtual inline void set_fe_var_data(const fe_var_type& fe_data)
    {
        libmesh_assert_msg(!_fe_var_data, "FE data already initialized.");
        _fe_var_data = &fe_data;
    }

    virtual inline uint_type n_dofs() const {

        libmesh_assert_msg(_fe_var_data, "FE data not initialized.");
        return 2*_fe_var_data->get_fe_shape_data().n_basis();
    }
    
    virtual inline void operator() (ContextType& c, vector_type& res, matrix_type* jac = nullptr) const {
        
        libmesh_assert_msg(_fe_var_data, "FE data not initialized.");
        libmesh_assert_msg(_property, "Section property not initialized");
        
        const typename fe_var_type::fe_shape_data_type
        &fe = _fe_var_data->get_fe_shape_data();
        
        typename EigenVector<var_scalar_type>::type
        epsilon = EigenVector<var_scalar_type>::type::Zero(3),
        stress  = EigenVector<var_scalar_type>::type::Zero(3),
        vec     = EigenVector<var_scalar_type>::type::Zero(2*fe.n_basis());
        
        typename EigenMatrix<var_scalar_type>::type
        mat     = EigenMatrix<var_scalar_type>::type::Zero(3, 3),
        mat1    = EigenMatrix<var_scalar_type>::type::Zero(3, 2*fe.n_basis()),
        mat2    = EigenMatrix<var_scalar_type>::type::Zero(2*fe.n_basis(), 2*fe.n_basis());

        MAST::FEMOperatorMatrixBase<nodal_scalar_type>
        Bxmat;
        Bxmat.reinit(3, 2, fe.n_basis()); // three stress-strain components

        
        for (uint_type i=0; i<fe.n_q_points(); i++) {
            
            c.qp = i;
            
            _property->value(c, mat);
            MAST::Element2D::linear_continuum_strain<nodal_scalar_type, var_scalar_type, fe_var_type>(*_fe_var_data, i, epsilon, Bxmat);
            stress = mat * epsilon;
            Bxmat.vector_mult_transpose(vec, stress);
            res += fe.detJxW(i) * vec;
            
            if (jac) {
                
                Bxmat.left_multiply(mat1, mat);
                Bxmat.right_multiply_transpose(mat2, mat1);
                (*jac) += fe.detJxW(i) * mat2;
            }
        }
    }
    
protected:
    
    
    const section_property_type       *_property;
    const fe_var_type                 *_fe_var_data;
};
}
}

#endif // _mast_strain_energy_compute_kernel_h_
