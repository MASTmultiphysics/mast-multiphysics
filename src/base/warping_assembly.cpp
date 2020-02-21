/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2019  Manav Bhatia
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

// MAST includes
#include "base/warping_assembly.h"
#include "base/system_initialization.h"
#include "base/physics_discipline_base.h"
#include "base/mesh_field_function.h"
#include "base/nonlinear_system.h"
#include "base/nonlinear_implicit_assembly_elem_operations.h"
#include "boundary_condition/point_load_condition.h"
#include "numerics/utility.h"
#include "mesh/geom_elem.h"
#include "mesh/fe_base.h"
#include "property_cards/element_property_card_2D.h"
#include "property_cards/material_property_card_base.h"

// libMesh includes
#include "libmesh/nonlinear_solver.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dof_map.h"
#include "libmesh/petsc_matrix.h"



MAST::WarpingAssembly::
WarpingAssembly():MAST::AssemblyBase(),
_post_assembly           (nullptr),
_res_l2_norm             (0.),
_first_iter_res_l2_norm  (-1.) {
    
}



MAST::WarpingAssembly::~WarpingAssembly() {
    
}



void
MAST::WarpingAssembly::
set_post_assembly_operation(MAST::WarpingAssembly::PostAssemblyOperation& post) {
    
    _post_assembly = &post;
}


void
MAST::WarpingAssembly::
make_matrix_symmetric(libMesh::SparseMatrix<Real>* J)
{
    libMesh::PetscMatrix<Real> Jt(_system->system().comm());
    J->close(); // Necessary to be able to do the transpose
    J->get_transpose(Jt);
    J->add(1.0, Jt);
    J->add(-0.5, *J);
}


const geometric_properties MAST::WarpingAssembly::calculate_geometric_properties() const
{
    geometric_properties gp;
    
    MAST::NonlinearSystem& nonlin_sys = _system->system();
        
    libMesh::MeshBase::const_element_iterator el = nonlin_sys.get_mesh().active_local_elements_begin();
        
    const libMesh::MeshBase::const_element_iterator end_el = nonlin_sys.get_mesh().active_local_elements_end();
    
    for ( ; el != end_el; ++el) 
    {
        const libMesh::Elem* elem = *el;
        
        uint n_nodes = elem->n_nodes();
        
        MAST::GeomElem geom_elem;
        geom_elem.init(*elem, *_system);
        
        /** Default number of quadrature points for second order triangular 
         * elements is 7.  We only need six for full integration of section 
         * properties. So, we set the number of extra quadrature points to -1 
         * for 6 quadrature points
         */
        const MAST::ElementPropertyCard2D& property =
            dynamic_cast<const MAST::ElementPropertyCard2D&>(_discipline->get_property_card(*elem));
            
        std::unique_ptr<MAST::FEBase> fe(geom_elem.init_fe(true, false, property.extra_quadrature_order(geom_elem)));
        
        const std::vector<Real>& JxW                    = fe->get_JxW();
        const std::vector<libMesh::Point>& qpoint       = fe->get_xyz();
        const std::vector<std::vector<Real>>& phi       = fe->get_phi();
        //const std::vector<std::vector<libMesh::RealVectorValue>>& dphi = fe->get_dphi();
        
        gp.A += elem->volume();
        
        // Get column vectors of x and y coordinates
//         RealVectorX xe =  RealVectorX::Zero(n_nodes);
//         RealVectorX ye =  RealVectorX::Zero(n_nodes);
//         for (uint i=0; i<n_nodes; i++) // Loop over nodes belonging to element
//         {
//             xe(i) = elem->point(i)(0);
//             ye(i) = elem->point(i)(1);
//         }
        
        // Loop over quadrature points
        for (unsigned int qp=0; qp<qpoint.size(); qp++) 
        {
            RealVectorX N = RealVectorX::Zero(n_nodes);
            
//             for (uint i=0; i<n_nodes; i++)
//             {
//                 N(i) = phi[i][qp];
//             }
//             Real Nye = N.dot(ye);
//             Real Nxe = N.dot(xe);
            
            Real Nye = qpoint[qp](1);
            Real Nxe = qpoint[qp](0);
            
            gp.Qx += JxW[qp] * Nye;
            gp.Qy += JxW[qp] * Nxe;
            
            gp.Ixx += JxW[qp] * Nye * Nye;
            gp.Iyy += JxW[qp] * Nxe * Nxe;
            gp.Ixy += JxW[qp] * Nye * Nxe;
        } // end loop over quadrature points
    } // end loop over elements

    gp.xc = gp.Qy/gp.A;
    gp.yc = gp.Qx/gp.A;

    gp.Ixxc = gp.Ixx - gp.Qx*gp.Qx/gp.A;
    gp.Iyyc = gp.Iyy - gp.Qy*gp.Qy/gp.A;
    gp.Ixyc = gp.Ixy - gp.Qx*gp.Qy/gp.A;
    
    gp.Ip = gp.Ixxc + gp.Iyyc; // TODO: double check this

    gp.rx = sqrt(gp.Ixx/gp.A);
    gp.ry = sqrt(gp.Iyy/gp.A);

    Real delta = sqrt( pow((gp.Ixxc - gp.Iyyc)/2.0, 2.0) + gp.Ixyc*gp.Ixyc );
    gp.I11 = (gp.Ixxc+gp.Iyyc)/2.0 + delta;
    gp.I22 = (gp.Ixxc+gp.Iyyc)/2.0 - delta;
    gp.phi_p = atan2(gp.Ixxc-gp.I11, gp.Ixyc);
    
    return gp;
}


const warping_properties MAST::WarpingAssembly::calculate_warping_properties(
    const libMesh::NumericVector<Real>& F_warp, 
    const libMesh::NumericVector<Real>& Omega, 
    const libMesh::NumericVector<Real>& Psi, 
    const libMesh::NumericVector<Real>& Phi, 
    const Real A, const Real Ixxc, const Real Iyyc, const Real Ixyc,
    const Real xc, const Real yc) const
{
    warping_properties wp;
    
    MAST::NonlinearSystem& nonlin_sys = _system->system();
    
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = nonlin_sys.get_dof_map();
        
    const uint w_var = _system->system().variable_number("warping_w");
    const uint lambda_var = _system->system().variable_number("warping_lambda");
    
    libMesh::MeshBase::const_element_iterator el = nonlin_sys.get_mesh().active_local_elements_begin();
        
    const libMesh::MeshBase::const_element_iterator end_el = nonlin_sys.get_mesh().active_local_elements_end();
    
    Real kappa_x  = 0.0;
    Real kappa_y  = 0.0;
    Real kappa_xy = 0.0;
    
    Real inv_kappa_x  = 0.0;
    Real inv_kappa_y  = 0.0;
    Real inv_kappa_xy = 0.0;
    
    Real nu, delta_s;
    
    for ( ; el != end_el; ++el) 
    {
        const libMesh::Elem* elem = *el;
        
        uint n_nodes = elem->n_nodes();
        
        MAST::GeomElem geom_elem;
        geom_elem.init(*elem, *_system);
        
        dof_map.dof_indices(elem, dof_indices, w_var);
        
        const MAST::ElementPropertyCard2D& property =
            dynamic_cast<const MAST::ElementPropertyCard2D&>(_discipline->get_property_card(*elem));
            
        const MAST::FieldFunction<Real>& nu_f = property.get_material().get<MAST::FieldFunction<Real>>("nu");
        
        /** Default number of quadrature points for second order triangular 
         * elements is 7.  We only need six for full integration of section 
         * properties. So, we set the number of extra quadrature points to -1 
         * for 6 quadrature points
         */
        std::unique_ptr<MAST::FEBase> fe(geom_elem.init_fe(true, false, property.extra_quadrature_order(geom_elem)));
        
        const std::vector<Real>& JxW                    = fe->get_JxW();
        const std::vector<libMesh::Point>& qpoint       = fe->get_xyz();
        const std::vector<std::vector<Real>>& phi       = fe->get_phi();
        const std::vector<std::vector<libMesh::RealVectorValue>>& dphi = fe->get_dphi();
                
        // Get column vectors of x and y coordinates
//         RealVectorX xe =  RealVectorX::Zero(n_nodes);
//         RealVectorX ye =  RealVectorX::Zero(n_nodes);
        RealVectorX we =  RealVectorX::Zero(n_nodes);
        RealVectorX psie =  RealVectorX::Zero(n_nodes);
        RealVectorX phie =  RealVectorX::Zero(n_nodes);
        for (uint i=0; i<n_nodes; i++) // Loop over nodes belonging to element
        {
//             xe(i) = elem->point(i)(0);
//             ye(i) = elem->point(i)(1);
            we(i) = Omega(dof_indices[i]);
            psie(i) = Psi(dof_indices[i]);
            phie(i) = Phi(dof_indices[i]);
        }
        
        // Loop over quadrature points
        for (unsigned int qp=0; qp<qpoint.size(); qp++) 
        {
            RealVectorX N = RealVectorX::Zero(n_nodes);
            RealMatrixX B = RealMatrixX::Zero(2, n_nodes);
            for (uint i=0; i<n_nodes; i++)
            {
                N(i) = phi[i][qp];
                B(0,i) = dphi[i][qp](0);
                B(1,i) = dphi[i][qp](1);
            }
            
            nu_f(qpoint[qp], 0.0, nu); // NOTE: Hard-coded time value of 0.0, does this need to be variable?
            
//             Real Nye = N.dot(ye);
//             Real Nxe = N.dot(xe);
            Real Nye = qpoint[qp](1);
            Real Nxe = qpoint[qp](0);
            Real Nwe = N.dot(we);
            
            wp.Qw += JxW[qp] * Nwe;
            wp.Iw += JxW[qp] * Nwe*Nwe;
            
            wp.Ixw += JxW[qp] * Nxe * Nwe;
            wp.Iyw += JxW[qp] * Nye * Nwe;
            
            delta_s = 2.0*(1.0+nu)*(Ixxc*Iyyc - Ixyc*Ixyc);
            
            wp.xs += JxW[qp] * nu/(2.0*delta_s) * (Iyyc*Nxe + Ixyc*Nye) * (Nxe*Nxe + Nye*Nye);
            wp.ys += JxW[qp] * nu/(2.0*delta_s) * (Ixxc*Nye + Ixyc*Nxe) * (Nxe*Nxe + Nye*Nye);
            
            Real r = Nxe*Nxe - Nye*Nye;
            Real q = 2.0 * Nxe * Nye;
            RealVectorX d = RealVectorX::Zero(2);
            d(0) = Ixxc*r - Ixyc*q;
            d(1) = Ixyc*r + Ixxc*q;
            RealVectorX h = RealVectorX::Zero(2);
            h(0) = -Ixyc*r + Iyyc*q;
            h(1) = -Iyyc*r - Ixyc*q;
            
            RealVectorX v = B*psie - nu/2.0 * d;
            RealVectorX w = B*phie - nu/2.0 * h;
            inv_kappa_x  += JxW[qp] * v.transpose().dot(v) / (delta_s*delta_s);
            inv_kappa_y  += JxW[qp] * w.transpose().dot(w) / (delta_s*delta_s);
            inv_kappa_xy += JxW[qp] * v.transpose().dot(w) / (delta_s*delta_s);
        } // end loop over quadrature points
        dof_indices.clear();
    } // end loop over elements
    
    // TODO: This isn't quite correct below because it doesn't account for a 
    // changing value of nu like is done above when looping through elements
    // and quadrature points.
    wp.xs -= F_warp.dot(Phi)/delta_s;
    wp.ys += F_warp.dot(Psi)/delta_s;
//     std::unique_ptr<libMesh::NumericVector<Real>> scaled_F_warp = F_warp.clone();
//     uint i=0;
//     for (const auto node_ptr : nonlin_sys.get_mesh().node_ptr_range())
//     {
//         dof_map.dof_indices(node_ptr, dof_indices, w_var);
//         nu_f(*node_ptr, 0.0, nu);
//         delta_s = 2.0*(1.0+nu)*(Ixxc*Iyyc - Ixyc*Ixyc);
//         scaled_F_warp->operator()(i) /= delta_s;
//         dof_indices.clear();
//     }
//     wp.xs -= scaled_F_warp->dot(Phi);
//     wp.ys += scaled_F_warp->dot(Psi);
    
    Real denom = Ixxc*Iyyc - Ixyc*Ixyc;
    wp.xs_t = (Ixyc * wp.Ixw - Iyyc * wp.Iyw) / denom;
    wp.ys_t = (Ixxc * wp.Ixw - Ixyc * wp.Iyw) / denom;
    wp.gamma = wp.Iw - wp.Qw * wp.Qw / A - wp.ys * wp.Ixw + wp.xs * wp.Iyw;
    
    wp.kappa_x = 1.0/(A*inv_kappa_x);
    wp.kappa_y = 1.0/(A*inv_kappa_y);
    wp.kappa_xy = 1.0/(A*inv_kappa_xy);
    
    wp.J = Ixxc + Iyyc - Omega.dot(F_warp);
    
    // Shift the shear centers to account for the shifting of the mesh's
    // global coordinates to the sections geometric centroid
    wp.xs_t += xc;
    wp.ys_t += yc;
    wp.xs += xc;
    wp.ys += yc;
    
    return wp;
}


void 
MAST::WarpingAssembly::
get_loads(libMesh::NumericVector<Real>& F_warp, 
          libMesh::NumericVector<Real>& F_shearx,
          libMesh::NumericVector<Real>& F_sheary,
          const Real xc, const Real yc,
          const Real Ixxc, const Real Iyyc, const Real Ixyc)
{
    MAST::NonlinearSystem& nonlin_sys = _system->system();
    
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = nonlin_sys.get_dof_map();
    
    libMesh::MeshBase::const_element_iterator el = 
        nonlin_sys.get_mesh().active_local_elements_begin();
        
    const libMesh::MeshBase::const_element_iterator end_el = 
        nonlin_sys.get_mesh().active_local_elements_end();
    
    const uint w_var = _system->system().variable_number("warping_w");
    const uint lambda_var = _system->system().variable_number("warping_lambda");
    
    for ( ; el != end_el; ++el) 
    {
        const libMesh::Elem* elem = *el;
        
        MAST::GeomElem geom_elem;
        geom_elem.init(*elem, *_system);
        
        const MAST::ElementPropertyCard2D& property =
            dynamic_cast<const MAST::ElementPropertyCard2D&>(_discipline->get_property_card(*elem));
            
        const MAST::FieldFunction<Real>& nu_f = property.get_material().get<MAST::FieldFunction<Real>>("nu");
            
        dof_map.dof_indices(elem, dof_indices, w_var);
        
        /** Default number of quadrature points for second order triangular 
         * elements is 7.  We only need six for full integration of section 
         * properties. So, we set the number of extra quadrature points to -1 
         * for 6 quadrature points
         */
        std::unique_ptr<MAST::FEBase> fe(geom_elem.init_fe(true, false, property.extra_quadrature_order(geom_elem)));
        
        const std::vector<Real>& JxW                    = fe->get_JxW();
        const std::vector<libMesh::Point>& qpoint       = fe->get_xyz();
        const std::vector<std::vector<Real>>& phi       = fe->get_phi();
        const std::vector<std::vector<libMesh::RealVectorValue>>& dphi = fe->get_dphi();
                
        // Get column vectors of x and y coordinates of element nodes, shifted
        // so that structures centriod is at (0,0).
        uint n_nodes = elem->n_nodes();
//         RealVectorX xe =  RealVectorX::Zero(n_nodes);
//         RealVectorX ye =  RealVectorX::Zero(n_nodes);
//         for (uint i=0; i<n_nodes; i++) // Loop over nodes belonging to element
//         {
//             xe(i) = elem->node_ref(i)(0);// - xc;
//             ye(i) = elem->node_ref(i)(1);// - yc;
//         }
        
        RealVectorX f_warp_e = RealVectorX::Zero(n_nodes);
        RealVectorX f_shear_x_e = RealVectorX::Zero(n_nodes);
        RealVectorX f_shear_y_e = RealVectorX::Zero(n_nodes);
        
        // Loop over quadrature points
        for (unsigned int qp=0; qp<qpoint.size(); qp++) 
        {
            RealVectorX N = RealVectorX::Zero(n_nodes);
            RealMatrixX B = RealMatrixX::Zero(2, n_nodes);
            for (uint i=0; i<n_nodes; i++)
            {
                N(i) = phi[i][qp];
                B(0,i) = dphi[i][qp](0);
                B(1,i) = dphi[i][qp](1);
            }
            
//             Real Nye = N.dot(ye);
//             Real Nxe = N.dot(xe);
            Real Nye = qpoint[qp](1);
            Real Nxe = qpoint[qp](0);
            
            RealVectorX v = RealVectorX::Zero(2);
            v(0) = Nye;
            v(1) = -Nxe;
            
            f_warp_e += JxW[qp] * B.transpose() * v;
            
            Real r = Nxe*Nxe - Nye*Nye;
            Real q = 2.0 * Nxe * Nye;
            RealVectorX d = RealVectorX::Zero(2);
            d(0) = Ixxc*r - Ixyc*q;
            d(1) = Ixyc*r + Ixxc*q;
            RealVectorX h = RealVectorX::Zero(2);
            h(0) = -Ixyc*r + Iyyc*q;
            h(1) = -Iyyc*r - Ixyc*q;
            
            Real nu;
            nu_f(qpoint[qp], 0.0, nu); // NOTE: Hard-coded time of 0.0, will this ever need to be variable?
            
            f_shear_x_e += JxW[qp] * (  nu/2. * B.transpose() * d + 2.0*(1.0+nu) * N * (Ixxc*Nxe - Ixyc*Nye) );
            f_shear_y_e += JxW[qp] * (  nu/2. * B.transpose() * h + 2.0*(1.0+nu) * N * (Iyyc*Nye - Ixyc*Nxe) );
        } // end loop over quadrature points
        
        DenseRealVector v;
        
        MAST::copy(v, -f_warp_e);
        F_warp.add_vector(v, dof_indices);
        
        MAST::copy(v, f_shear_x_e);
        F_shearx.add_vector(v, dof_indices);
        
        MAST::copy(v, f_shear_y_e);
        F_sheary.add_vector(v, dof_indices);
        
        dof_indices.clear();
    } // end loop over elements
    
    F_warp.close();
    F_shearx.close();
    F_sheary.close();
} // end get_loads() function



// /**
//  *  Need w (omega) to calculate:
//  *      Jt, Ixw, Iyw, xs_t, ys_t, Qw, Iw, Gamma_t
//  * 
//  *  Need Psi and Theta to calculate:
//  *      xs, ys, Gamma, kappa_x, kappa_y, kappa_xy
//  */

void
MAST::WarpingAssembly::
residual_and_jacobian (const libMesh::NumericVector<Real>& X,
                       libMesh::NumericVector<Real>* R,
                       libMesh::SparseMatrix<Real>*  J,
                       libMesh::NonlinearImplicitSystem& S) {
    
    libmesh_assert(_system);
    libmesh_assert(_discipline);
    libmesh_assert(_elem_ops);

    MAST::NonlinearSystem& nonlin_sys = _system->system();
    
    // make sure that the system for which this object was created,
    // and the system passed through the function call are the same
    libmesh_assert_equal_to(&S, &(nonlin_sys));
    
    if (R) R->zero();
    if (J) J->zero();
        
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX vec, sol;
    RealMatrixX mat;
    
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = _system->system().get_dof_map();
    
    
    std::unique_ptr<libMesh::NumericVector<Real> > localized_solution;
    localized_solution.reset(build_localized_vector(nonlin_sys,
                                                     X).release());
    
    
    // if a solution function is attached, initialize it
    if (_sol_function)
        _sol_function->init( X);
    
    
    libMesh::MeshBase::const_element_iterator       el     =
    nonlin_sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    nonlin_sys.get_mesh().active_local_elements_end();
    
    MAST::NonlinearImplicitAssemblyElemOperations&
    ops = dynamic_cast<MAST::NonlinearImplicitAssemblyElemOperations&>(*_elem_ops);
    
    const uint w_var = _system->system().variable_number("warping_w");
    const uint lambda_var = _system->system().variable_number("warping_lambda");

    for ( ; el != end_el; ++el) 
    {
        const libMesh::Elem* elem = *el;
        
        dof_map.dof_indices (elem, dof_indices, w_var);
        
        MAST::GeomElem geom_elem;
        ops.set_elem_data(elem->dim(), *elem, geom_elem);
        geom_elem.init(*elem, *_system);
        
        ops.init(geom_elem);

        // get the solution
        unsigned int ndofs = (unsigned int)dof_indices.size();
        sol.setZero(ndofs);
        vec.setZero(ndofs);
        mat.setZero(ndofs, ndofs);
        
        for (unsigned int i=0; i<dof_indices.size(); i++)
            sol(i) = (*localized_solution)(dof_indices[i]);
        
        ops.set_elem_solution(sol);
        
        
//        if (_sol_function)
//            physics_elem->attach_active_solution_function(*_sol_function);
        
        //_check_element_numerical_jacobian(*physics_elem, sol);
        
        // perform the element level calculations
        ops.elem_calculations(J!=nullptr?true:false,
                                              vec, mat);
        
//        physics_elem->detach_active_solution_function();

        ops.clear_elem();
        
        // copy to the libMesh matrix for further processing
        DenseRealVector v;
        DenseRealMatrix m;
        if (R)
            MAST::copy(v, vec);
        if (J)
            MAST::copy(m, mat);

        // constrain the quantities to account for hanging dofs,
        // Dirichlet constraints, etc.
        if (R && J)
            dof_map.constrain_element_matrix_and_vector(m, v, dof_indices);
        else if (R)
            dof_map.constrain_element_vector(v, dof_indices);
        else
            dof_map.constrain_element_matrix(m, dof_indices);
        
        // add to the global matrices
        if (R) R->add_vector(v, dof_indices);
        if (J) J->add_matrix(m, dof_indices);
        dof_indices.clear();
    }
    
    if (J)
    {
        uint n = J->n();
        for (uint i=0; i<n-1; i++)
        {
            J->add(n-1,i,1.0);
            J->add(i,n-1,1.0);
        }
        //J->add(J->m()-1, J->n()-1, 1.5e-08);
    }
    
    // call the post assembly object, if provided by user
    if (_post_assembly)
        _post_assembly->post_assembly(X, R, J, S);
    

    // if a solution function is attached, clear it
    if (_sol_function)
        _sol_function->clear();
    
    if (R) {
        
        R->close();
        _res_l2_norm = R->l2_norm();
        if (_first_iter_res_l2_norm < 0.)
            _first_iter_res_l2_norm = _res_l2_norm;
    }
    
    if ((J) and (_force_jacobian_symmetry)) 
    {
        make_matrix_symmetric(J);
    }
    
    if (J && close_matrix) J->close();
}




void
MAST::WarpingAssembly::
linearized_jacobian_solution_product (const libMesh::NumericVector<Real>& X,
                                      const libMesh::NumericVector<Real>& dX,
                                      libMesh::NumericVector<Real>& JdX,
                                      libMesh::NonlinearImplicitSystem& S) 
{
    libmesh_error_msg("linearized_jacobian_solution_product not implemented in warping_assembly.cpp");
}



void
MAST::WarpingAssembly::
second_derivative_dot_solution_assembly (const libMesh::NumericVector<Real>& X,
                                         const libMesh::NumericVector<Real>& dX,
                                         libMesh::SparseMatrix<Real>& d_JdX_dX,
                                         libMesh::NonlinearImplicitSystem& S) 
{
    libmesh_error_msg("second_derivative_dot_solution_assembly not implemented in warping_assembly.cpp");
}





bool
MAST::WarpingAssembly::
sensitivity_assemble (const MAST::FunctionBase& f,
                      libMesh::NumericVector<Real>& sensitivity_rhs) 
{
    libmesh_error_msg("sensitivity_assemble not implemented in warping_assembly.cpp");
    
    libmesh_assert(_system);
    libmesh_assert(_discipline);
    libmesh_assert(_elem_ops);

    MAST::NonlinearSystem& nonlin_sys = _system->system();
    
    sensitivity_rhs.zero();
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX vec, sol;
    
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = nonlin_sys.get_dof_map();
    
    
    std::unique_ptr<libMesh::NumericVector<Real> > localized_solution;
    localized_solution.reset(build_localized_vector(nonlin_sys,
                                                     *nonlin_sys.solution).release());
    
    // if a solution function is attached, initialize it
    if (_sol_function)
        _sol_function->init( *nonlin_sys.solution);
    
    libMesh::MeshBase::const_element_iterator       el     =
    nonlin_sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    nonlin_sys.get_mesh().active_local_elements_end();
    
    MAST::NonlinearImplicitAssemblyElemOperations&
    ops = dynamic_cast<MAST::NonlinearImplicitAssemblyElemOperations&>(*_elem_ops);

    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem = *el;
        
        // no sensitivity computation assembly is neeed in these cases
        if (_param_dependence &&
            // if object is specified and elem does not depend on it
            !_param_dependence->if_elem_depends_on_parameter(*elem, f))
            continue;

        dof_map.dof_indices (elem, dof_indices);
        
        MAST::GeomElem geom_elem;
        ops.set_elem_data(elem->dim(), *elem, geom_elem);
        geom_elem.init(*elem, *_system);
        
        ops.init(geom_elem);

        // get the solution
        unsigned int ndofs = (unsigned int)dof_indices.size();
        sol.setZero(ndofs);
        vec.setZero(ndofs);
        
        for (unsigned int i=0; i<dof_indices.size(); i++)
            sol(i) = (*localized_solution)(dof_indices[i]);
        
        ops.set_elem_solution(sol);
        
//        if (_sol_function)
//            physics_elem->attach_active_solution_function(*_sol_function);
        
        ops.elem_sensitivity_calculations(f, vec);
        
//        physics_elem->detach_active_solution_function();
        ops.clear_elem();

        // copy to the libMesh matrix for further processing
        DenseRealVector v;
        MAST::copy(v, vec);

        // constrain the quantities to account for hanging dofs,
        // Dirichlet constraints, etc.
        dof_map.constrain_element_vector(v, dof_indices);
        
        // add to the global matrices
        sensitivity_rhs.add_vector(v, dof_indices);
        dof_indices.clear();
    }
    
    // add the point loads if any in the discipline
    if (_discipline->point_loads().size()) {
        
        const MAST::PointLoadSetType&
        loads = _discipline->point_loads();
        
        vec = RealVectorX::Zero(_system->n_vars());
        
        MAST::PointLoadSetType::const_iterator
        it    = loads.begin(),
        end   = loads.end();
        
        const libMesh::dof_id_type
        first_dof  = dof_map.first_dof(nonlin_sys.comm().rank()),
        end_dof	   = dof_map.end_dof(nonlin_sys.comm().rank());
        
        for ( ; it != end; it++) {
            
            // get the point load function
            const MAST::FieldFunction<RealVectorX>
            &func = (*it)->get<MAST::FieldFunction<RealVectorX>>("load");
            
            // get the nodes on which this object defines the load
            const std::set<const libMesh::Node*>
            nodes = (*it)->get_nodes();
            
            std::set<const libMesh::Node*>::const_iterator
            n_it    = nodes.begin(),
            n_end   = nodes.end();
            
            for (; n_it != n_end; n_it++) {
                
                // load at the node
                vec.setZero();
                func.derivative(f, **n_it, nonlin_sys.time, vec);
                vec *= -1.;
                
                dof_map.dof_indices(*n_it, dof_indices);
                
                libmesh_assert_equal_to(dof_indices.size(), vec.rows());

                // zero the components of the vector if they do not
                // belong to this processor
                for (unsigned int i=0; i<dof_indices.size(); i++)
                    if (dof_indices[i] <   first_dof  ||
                        dof_indices[i] >=  end_dof)
                        vec(i) = 0.;

                DenseRealVector v;
                MAST::copy(v, vec);
                
                dof_map.constrain_element_vector(v, dof_indices);
                sensitivity_rhs.add_vector(v, dof_indices);
                dof_indices.clear();
            }
        }
    }

    // if a solution function is attached, initialize it
    if (_sol_function)
        _sol_function->clear();
    
    sensitivity_rhs.close();
    
    return true;
}

void MAST::WarpingAssembly::set_force_jacobian_symmetry(bool tf)
{
    _force_jacobian_symmetry = tf;
}

const bool MAST::WarpingAssembly::get_force_jacobian_symmetry() const
{
    return _force_jacobian_symmetry;
}
