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
#include "fluid/frequency_domain_linearized_conservative_fluid_elem.h"
#include "fluid/primitive_fluid_solution.h"
#include "fluid/small_disturbance_primitive_fluid_solution.h"
#include "fluid/flight_condition.h"
#include "base/boundary_condition_base.h"
#include "base/system_initialization.h"
#include "aeroelasticity/frequency_function.h"
#include "elasticity/normal_rotation_function_base.h"
#include "base/nonlinear_system.h"
#include "mesh/fe_base.h"
#include "mesh/geom_elem.h"
#include "base/assembly_base.h"



MAST::FrequencyDomainLinearizedConservativeFluidElem::
FrequencyDomainLinearizedConservativeFluidElem(MAST::SystemInitialization& sys,
                                               MAST::AssemblyBase& assembly,
                                               const MAST::GeomElem& elem,
                                               const MAST::FlightCondition& f):
MAST::ConservativeFluidElementBase(sys, assembly, elem, f),
freq(nullptr) {
    
    
}





MAST::FrequencyDomainLinearizedConservativeFluidElem::
~FrequencyDomainLinearizedConservativeFluidElem() {
    
    
}




bool
MAST::FrequencyDomainLinearizedConservativeFluidElem::
internal_residual (bool request_jacobian,
                   ComplexVectorX& f,
                   ComplexMatrixX& jac) {
    
    // first get the internal residual and Jacobian from the
    // conservative elems, and then add to it an extra dissipation
    // to control solution discontinuities emanating from the boundary.
    
    //const std::vector<Real>& JxW                  = _fe->get_JxW();
    //const std::vector<std::vector<Real> >& phi    = _fe->get_phi();
    const unsigned int
    n2     = f.size();
    //nphi   = _fe->n_shape_functions();
    
    RealMatrixX
    f_jac_x         = RealMatrixX::Zero(   n2,    n2),
    fm_jac_xdot     = RealMatrixX::Zero(   n2,    n2);
    /*mat1_n1n1       = RealMatrixX::Zero(   n1,    n1),
    mat2_n1n1       = RealMatrixX::Zero(   n1,    n1),
    mat3_n1n2       = RealMatrixX::Zero(   n1,    n2),
    mat4_n2n2       = RealMatrixX::Zero(   n2,    n2),
    AiBi_adv        = RealMatrixX::Zero(   n1,    n2),
    A_sens          = RealMatrixX::Zero(   n1,    n2),
    LS              = RealMatrixX::Zero(   n1,    n2),
    LS_sens         = RealMatrixX::Zero(   n2,    n2);*/
    
    
    ComplexMatrixX
    local_jac       = ComplexMatrixX::Zero(   n2,    n2);
    
    
    RealVectorX
    //vec1_n1  = RealVectorX::Zero(n1),
    //vec2_n1  = RealVectorX::Zero(n1),
    local_f    = RealVectorX::Zero(n2);
    //dc       = RealVectorX::Zero(dim);
    
    const Complex
    iota(0., 1.);

    Real
    omega   = 0.,
    b_V     = 0.;
    
    (*freq)(omega);
    freq->nondimensionalizing_factor(b_V);
    
    // df/dx. We always need the Jacobian, since it is used to calculate
    // the residual
    MAST::ConservativeFluidElementBase::internal_residual(true,
                                                          local_f,
                                                          f_jac_x);
    
    // dfm/dxdot. We always need the Jacobian, since it is used to calculate
    // the residual
    MAST::ConservativeFluidElementBase::velocity_residual(true,
                                                          local_f,
                                                          fm_jac_xdot,
                                                          f_jac_x);
    
    // now, combine the two to return the complex Jacobian
    
    local_jac = (f_jac_x.cast<Complex>() * b_V +  // multiply stiffness with nondim factor
                 iota * omega * fm_jac_xdot.cast<Complex>());
    
    if (request_jacobian)
        jac      +=  local_jac;
    
    f        +=  local_jac * _complex_sol;
    
    return request_jacobian;
}





bool
MAST::FrequencyDomainLinearizedConservativeFluidElem::
internal_residual_sensitivity (const MAST::FunctionBase& p,
                               bool request_jacobian,
                               ComplexVectorX& f,
                               ComplexMatrixX& jac) {
    
    // first get the internal residual and Jacobian from the
    // conservative elems, and then add to it an extra dissipation
    // to control solution discontinuities emanating from the boundary.
    
    //const std::vector<Real>& JxW                  = _fe->get_JxW();
    //const std::vector<std::vector<Real> >& phi    = _fe->get_phi();
    const unsigned int
    n2     = f.size();
    //nphi   = _fe->n_shape_functions();
    
    RealMatrixX
    f_jac_x         = RealMatrixX::Zero(   n2,    n2),
    fm_jac_xdot     = RealMatrixX::Zero(   n2,    n2);
    /*mat1_n1n1       = RealMatrixX::Zero(   n1,    n1),
     mat2_n1n1       = RealMatrixX::Zero(   n1,    n1),
     mat3_n1n2       = RealMatrixX::Zero(   n1,    n2),
     mat4_n2n2       = RealMatrixX::Zero(   n2,    n2),
     AiBi_adv        = RealMatrixX::Zero(   n1,    n2),
     A_sens          = RealMatrixX::Zero(   n1,    n2),
     LS              = RealMatrixX::Zero(   n1,    n2),
     LS_sens         = RealMatrixX::Zero(   n2,    n2);*/
    
    
    ComplexMatrixX
    local_jac       = ComplexMatrixX::Zero(   n2,    n2),
    local_jac_sens  = ComplexMatrixX::Zero(   n2,    n2);
    
    
    RealVectorX
    //vec1_n1  = RealVectorX::Zero(n1),
    //vec2_n1  = RealVectorX::Zero(n1),
    local_f    = RealVectorX::Zero(n2);
    //dc       = RealVectorX::Zero(dim);
    
    const Complex
    iota(0., 1.);
    
    Real
    omega   = 0.,
    domega  = 0.,
    b_V     = 0.;
    
    (*freq)(omega);
    freq->derivative(p, domega);
    freq->nondimensionalizing_factor(b_V);
    
    
    // df/dx. We always need the Jacobian, since it is used to calculate
    // the residual
    MAST::ConservativeFluidElementBase::internal_residual(true,
                                                          local_f,
                                                          f_jac_x);
    
    // dfm/dxdot. We always need the Jacobian, since it is used to calculate
    // the residual
    MAST::ConservativeFluidElementBase::velocity_residual(true,
                                                          local_f,
                                                          fm_jac_xdot,
                                                          f_jac_x);
    
    // now, combine the two to return the complex Jacobian
    
    local_jac      = (f_jac_x.cast<Complex>() * b_V +  // multiply stiffness with nondim factor
                      iota * omega * fm_jac_xdot.cast<Complex>());

    local_jac_sens = (iota * domega * fm_jac_xdot.cast<Complex>());

    
    if (request_jacobian)
        jac      +=  local_jac_sens;
    
    f        +=  local_jac_sens * _complex_sol + local_jac * _complex_sol_sens;
    
    return request_jacobian;
}




bool
MAST::FrequencyDomainLinearizedConservativeFluidElem::
side_external_residual (bool request_jacobian,
                        ComplexVectorX& f,
                        ComplexMatrixX& jac,
                        std::multimap<libMesh::boundary_id_type, MAST::BoundaryConditionBase*>& bc) {

    // the far-field and symmetry wall, which are assumed to be stationary will
    // only contribute to the real part of the complex Jacobian. Hence,
    // we will use the parent class's methods for these two. The
    // slip wall, which may be oscillating, is implemented for this element.
    
    const unsigned int
    dim     = _elem.dim(),
    n1      = dim+2,
    n2      = f.size();
    
    RealMatrixX
    f_jac_x = RealMatrixX::Zero(n2, n2);
    
    RealVectorX
    local_f = RealVectorX::Zero(n2);
    
    Real
    b_V     = 0.;
    freq->nondimensionalizing_factor(b_V);

    
    std::map<unsigned int, std::vector<MAST::BoundaryConditionBase*>> loads;
    _elem.external_side_loads_for_quadrature_elem(bc, loads);
    
    std::map<unsigned int, std::vector<MAST::BoundaryConditionBase*>>::const_iterator
    it   = loads.begin(),
    end  = loads.end();
    
    for ( ; it != end; it++) {
        
        std::vector<MAST::BoundaryConditionBase*>::const_iterator
        bc_it  = it->second.begin(),
        bc_end = it->second.end();
        
        for ( ; bc_it != bc_end; bc_it++) {
            
            // apply all the types of loading
            switch ((*bc_it)->type()) {
                case MAST::SYMMETRY_WALL: {
                    
                    f_jac_x.setZero();
                    local_f.setZero();
                    
                    // We always need the Jacobian, since it is used to
                    // calculate the residual
                    
                    MAST::ConservativeFluidElementBase::
                    symmetry_surface_residual(true,
                                              local_f,
                                              f_jac_x,
                                              it->first,
                                              **bc_it);
                    
                    // multiply jacobian with the nondimensionalizing
                    // factor (V/b for flutter analysis)
                    if (request_jacobian)
                        jac  += f_jac_x.cast<Complex>() * b_V;
                    f    += f_jac_x.cast<Complex>() * b_V * _complex_sol;
                }
                    break;
                    
                case MAST::SLIP_WALL: {
                    
                    // this calculates the Jacobian and residual contribution
                    // including the nondimensionalizing factor.
                    this->slip_wall_surface_residual(request_jacobian,
                                                     f,
                                                     jac,
                                                     it->first,
                                                     **bc_it);
                }
                    break;
                    
                case MAST::FAR_FIELD: {
                    
                    f_jac_x.setZero();
                    local_f.setZero();
                    
                    // We always need the Jacobian, since it is used to
                    // calculate the residual
                    
                    MAST::ConservativeFluidElementBase::
                    far_field_surface_residual(true,
                                               local_f,
                                               f_jac_x,
                                               it->first,
                                               **bc_it);
                    
                    // multiply jacobian with the nondimensionalizing
                    // factor (V/b for flutter analysis)
                    if (request_jacobian)
                        jac  += f_jac_x.cast<Complex>() * b_V;
                    
                    f    += f_jac_x.cast<Complex>() * b_V * _complex_sol;
                }
                    break;
                    
                case MAST::DIRICHLET:
                    // nothing to be done here
                    break;
                    
                default:
                    // not implemented yet
                    libmesh_error();
                    break;
            }
        }
    }
    
    
    return request_jacobian;
}






bool
MAST::FrequencyDomainLinearizedConservativeFluidElem::
side_external_residual_sensitivity (const MAST::FunctionBase& p,
                                    bool request_jacobian,
                                    ComplexVectorX& f,
                                    ComplexMatrixX& jac,
                                    std::multimap<libMesh::boundary_id_type, MAST::BoundaryConditionBase*>& bc) {
    
    // the far-field and symmetry wall, which are assumed to be stationary will
    // only contribute to the real part of the complex Jacobian. Hence,
    // we will use the parent class's methods for these two. The
    // slip wall, which may be oscillating, is implemented for this element.
    
    const unsigned int
    dim     = _elem.dim(),
    n1      = dim+2,
    n2      = f.size();
    
    RealMatrixX
    f_jac_x = RealMatrixX::Zero(n2, n2);
    
    RealVectorX
    local_f = RealVectorX::Zero(n2);
    
    Real
    b_V     = 0.;
    freq->nondimensionalizing_factor(b_V);
    
    std::map<unsigned int, std::vector<MAST::BoundaryConditionBase*>> loads;
    _elem.external_side_loads_for_quadrature_elem(bc, loads);
    
    std::map<unsigned int, std::vector<MAST::BoundaryConditionBase*>>::const_iterator
    it   = loads.begin(),
    end  = loads.end();
    
    for ( ; it != end; it++) {
        
        std::vector<MAST::BoundaryConditionBase*>::const_iterator
        bc_it  = it->second.begin(),
        bc_end = it->second.end();
        
        for ( ; bc_it != bc_end; bc_it++) {
            
            // apply all the types of loading
            switch ((*bc_it)->type()) {
                case MAST::SYMMETRY_WALL: {
                    
                    f_jac_x.setZero();
                    local_f.setZero();
                    
                    // We always need the Jacobian, since it is used to
                    // calculate the residual
                    
                    MAST::ConservativeFluidElementBase::
                    symmetry_surface_residual(true,
                                              local_f,
                                              f_jac_x,
                                              it->first,
                                              **bc_it);
                    
                    // multiply jacobian with the nondimensionalizing
                    // factor (V/b for flutter analysis)
                    //if (request_jacobian)
                    //    jac  += f_jac_x.cast<Complex>() * b_V;
                    f    += f_jac_x.cast<Complex>() * b_V * _complex_sol_sens;
                }
                    break;
                    
                case MAST::SLIP_WALL: {
                    
                    // this calculates the Jacobian and residual contribution
                    // including the nondimensionalizing factor.
                    this->slip_wall_surface_residual_sensitivity(p, request_jacobian,
                                                                 f,
                                                                 jac,
                                                                 it->first,
                                                                 **bc_it);
                }
                    break;
                    
                case MAST::FAR_FIELD: {
                    
                    f_jac_x.setZero();
                    local_f.setZero();
                    
                    // We always need the Jacobian, since it is used to
                    // calculate the residual
                    
                    MAST::ConservativeFluidElementBase::
                    far_field_surface_residual(true,
                                               local_f,
                                               f_jac_x,
                                               it->first,
                                               **bc_it);
                    
                    // multiply jacobian with the nondimensionalizing
                    // factor (V/b for flutter analysis)
                    //if (request_jacobian)
                    //    jac  += f_jac_x.cast<Complex>() * b_V;
                    
                    f    += f_jac_x.cast<Complex>() * b_V * _complex_sol_sens;
                }
                    break;
                    
                case MAST::DIRICHLET:
                    // nothing to be done here
                    break;
                    
                default:
                    // not implemented yet
                    libmesh_error();
                    break;
            }
        }
    }
    
    return request_jacobian;
}






bool
MAST::FrequencyDomainLinearizedConservativeFluidElem::
slip_wall_surface_residual(bool request_jacobian,
                           ComplexVectorX& f,
                           ComplexMatrixX& jac,
                           const unsigned int s,
                           MAST::BoundaryConditionBase& bc) {
    
    // inviscid boundary condition without any diffusive component
    // conditions enforced are
    // vi ni = wi_dot (ni + dni) - ui dni   (moving slip wall with deformation)
    // tau_ij nj = 0   (because velocity gradient at wall = 0)
    // qi ni = 0       (since heat flux occurs only on no-slip wall and far-field bc)
    
    // prepare the side finite element
    std::unique_ptr<MAST::FEBase> fe(_elem.init_side_fe(s, false, false));

    const std::vector<Real> &JxW                 = fe->get_JxW();
    const std::vector<libMesh::Point>& normals   = fe->get_normals_for_reference_coordinate();
    const std::vector<libMesh::Point>& qpoint    = fe->get_xyz();
    
    const unsigned int
    dim        = _elem.dim(),
    n1         = dim+2,
    n2         = fe->n_shape_functions()*n1;
    
    RealVectorX
    vec1_n1    = RealVectorX::Zero(n1),
    uvec       = RealVectorX::Zero(3),
    dwdot_i    = RealVectorX::Zero(3),
    ni         = RealVectorX::Zero(3),
    dni        = RealVectorX::Zero(3),
    tmp        = RealVectorX::Zero(6);

    ComplexVectorX
    Dw_i       = ComplexVectorX::Zero(3),
    Dni        = ComplexVectorX::Zero(3),
    Duvec      = ComplexVectorX::Zero(3),
    vec2_n1    = ComplexVectorX::Zero(n1),
    vec2_n2    = ComplexVectorX::Zero(n2),
    flux       = ComplexVectorX::Zero(n1),
    tmp_c      = ComplexVectorX::Zero(6);

    RealMatrixX
    mat1_n1n1  = RealMatrixX::Zero( n1, n1),
    mat2_n1n2  = RealMatrixX::Zero( n1, n2),
    mat3_n2n2  = RealMatrixX::Zero( n2, n2);
    
    libMesh::Point pt;
    MAST::FEMOperatorMatrix Bmat;
    
    // create objects to calculate the primitive solution, flux, and Jacobian
    MAST::PrimitiveSolution                            primitive_sol;
    MAST::SmallPerturbationPrimitiveSolution<Complex>  sd_primitive_sol;
    
    Complex
    iota            (0., 1.),
    Dvi_ni_freq_indep   = 0.,
    Dvi_ni_freq_dep     = 0.;
    
    Real
    omega               = 0.,
    b_V                 = 0.,
    ui_ni_steady        = 0.;

    
    // get the surface motion object from the boundary condition object
    MAST::FieldFunction<RealVectorX>
    *vel   = nullptr;
    MAST::NormalRotationFunctionBase<RealVectorX>
    *n_rot = nullptr;

    MAST::FieldFunction<ComplexVectorX>
    *displ_perturb = nullptr;
    MAST::NormalRotationFunctionBase<ComplexVectorX>
    *n_rot_perturb = nullptr;

    if (bc.contains("velocity"))
        vel = &bc.get<MAST::FieldFunction<RealVectorX> >("velocity");
    
    if (bc.contains("normal_rotation")) {

        MAST::FieldFunction<RealVectorX>&
        tmp = bc.get<MAST::FieldFunction<RealVectorX> >("normal_rotation");
        n_rot = dynamic_cast<MAST::NormalRotationFunctionBase<RealVectorX>*>(&tmp);
    }

    if (bc.contains("frequency_domain_displacement")) {
        
        displ_perturb = &bc.get<MAST::FieldFunction<ComplexVectorX> >("frequency_domain_displacement");
        // if displ is provided then n_rot must also be provided
        libmesh_assert(bc.contains("frequency_domain_normal_rotation"));
        
        MAST::FieldFunction<ComplexVectorX>&
        tmp = bc.get<MAST::FieldFunction<ComplexVectorX> >("frequency_domain_normal_rotation");
        n_rot_perturb = dynamic_cast<MAST::NormalRotationFunctionBase<ComplexVectorX>*>(&tmp);
    }

    

    
    (*freq)(omega);
    freq->nondimensionalizing_factor(b_V);
    
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        // initialize the Bmat operator for this term
        _initialize_fem_interpolation_operator(qp, dim, *fe, Bmat);
        Bmat.right_multiply(vec1_n1,         _sol);  // conservative sol
        Bmat.right_multiply(vec2_n1, _complex_sol);  // perturbation sol
        
        // initialize the primitive solution
        primitive_sol.zero();
        primitive_sol.init(dim,
                           vec1_n1,
                           flight_condition->gas_property.cp,
                           flight_condition->gas_property.cv,
                           if_viscous());

        // initialize the small-disturbance primitive sol
        sd_primitive_sol.zero();
        sd_primitive_sol.init(primitive_sol, vec2_n1, if_viscous());
        
        ////////////////////////////////////////////////////////////
        //   Calculation of the surface velocity term. For a
        //   steady-state system,
        //        vi (ni + dni)  =  wdot_i (ni + dni)
        //  or    vi ni = wdot_i (ni + dni) - vi dni
        //
        //    The first order perturbation, D(.), of this is
        //        (vi + Dvi) (ni + dni + Dni) =
        //          (wdot_i + Dwdot_i) (ni + dni + Dni)
        //  or    Dvi ni = Dwdot_i (ni + dni + Dni) +
        //                  wdot_i (ni + dni + Dni) -
        //                  vi (ni + dni + Dni) -
        //                  Dvi (dni + Dni)
        //     Neglecting HOT, this becomes
        //        Dvi ni = Dwdot_i (ni + dni) +
        //                  wdot_i (ni + dni + Dni) -
        //                  vi (ni + dni + Dni) -
        //                  Dvi dni
        //               = Dwdot_i (ni + dni) +
        //                 wdot_i Dni - vi Dni - Dvi dni
        //
        ////////////////////////////////////////////////////////////
        
        // copy the surface normal
        for (unsigned int i_dim=0; i_dim<dim; i_dim++)
            ni(i_dim) = normals[qp](i_dim);
        
        // now check if the surface deformation is defined and
        // needs to be applied through transpiration boundary
        // condition
        primitive_sol.get_uvec(uvec);
        sd_primitive_sol.get_duvec(Duvec);
        
        flux.setZero();
        
        //////////////////////////////////////////////////////////////
        // contribution from the base-flow boundary condition
        //////////////////////////////////////////////////////////////
        if (vel) {
            
            (*vel)(qpoint[qp], _time, tmp);
            dwdot_i = tmp.topRows(3);
        }
        
        if (n_rot)
            (*n_rot)(qpoint[qp], normals[qp], _time, dni);

        ui_ni_steady  =  dwdot_i.dot(ni+dni) - uvec.dot(dni);
        
        flux         += ui_ni_steady * b_V * vec2_n1;               // vi_ni  dcons_flux
        flux(n1-1)   += ui_ni_steady * b_V * sd_primitive_sol.dp;   // vi_ni {0,0,0,0,Dp}
        
        if (displ_perturb) {
            
            //////////////////////////////////////////////////////////////
            // contribution from the small-disturbance boundary condition
            //////////////////////////////////////////////////////////////
            (*displ_perturb)(qpoint[qp], _time, tmp_c);
            Dw_i = tmp_c.topRows(3);
            (*n_rot_perturb)(qpoint[qp], normals[qp], _time, Dni);
        }
    
        Dvi_ni_freq_dep   = Dw_i.dot(ni+dni) * iota * omega;
        Dvi_ni_freq_indep = (dwdot_i.cast<Complex>().dot(Dni) -
                             uvec.cast<Complex>().dot(Dni) -
                             Duvec.dot(dni));
        
        
        flux             += (Dvi_ni_freq_indep * b_V +
                             Dvi_ni_freq_dep) * vec1_n1;              // Dvi_ni cons_flux
        flux(n1-1)       += (Dvi_ni_freq_indep * b_V +
                             Dvi_ni_freq_dep) * primitive_sol.p;      // Dvi_ni {0,0,0,0,p}
        
        // ni Dp  (only for the momentun eq)
        flux.segment(1, dim) += ((sd_primitive_sol.dp * b_V) *
                                 ni.segment(0,dim).cast<Complex>());
        
        Bmat.vector_mult_transpose(vec2_n2, flux);
        f += JxW[qp] * vec2_n2;
        
        if ( request_jacobian ) {
            
            // the Jacobian contribution comes only from the dp term in the
            // residual
            this->calculate_advection_flux_jacobian_for_moving_solid_wall_boundary
            (primitive_sol,
             ui_ni_steady,
             normals[qp],
             dni,
             mat1_n1n1);
            
            Bmat.left_multiply(mat2_n1n2, mat1_n1n1);
            Bmat.right_multiply_transpose(mat3_n2n2, mat2_n1n2);
            jac += (JxW[qp] * b_V) * mat3_n2n2.cast<Complex>() ;
        }
    }
    
    return request_jacobian;
}






bool
MAST::FrequencyDomainLinearizedConservativeFluidElem::
slip_wall_surface_residual_sensitivity(const MAST::FunctionBase& p,
                                       bool request_jacobian,
                                       ComplexVectorX& f,
                                       ComplexMatrixX& jac,
                                       const unsigned int s,
                                       MAST::BoundaryConditionBase& bc) {
    
    // inviscid boundary condition without any diffusive component
    // conditions enforced are
    // vi ni = wi_dot (ni + dni) - ui dni   (moving slip wall with deformation)
    // tau_ij nj = 0   (because velocity gradient at wall = 0)
    // qi ni = 0       (since heat flux occurs only on no-slip wall and far-field bc)
    
    // prepare the side finite element
    std::unique_ptr<MAST::FEBase> fe(_elem.init_side_fe(s, false, false));

    const std::vector<Real> &JxW                 = fe->get_JxW();
    const std::vector<libMesh::Point>& normals   = fe->get_normals_for_reference_coordinate();
    const std::vector<libMesh::Point>& qpoint    = fe->get_xyz();
    
    const unsigned int
    dim        = _elem.dim(),
    n1         = dim+2,
    n2         = fe->n_shape_functions()*n1;
    
    RealVectorX
    vec1_n1    = RealVectorX::Zero(n1),
    uvec       = RealVectorX::Zero(3),
    ni         = RealVectorX::Zero(3),
    dni        = RealVectorX::Zero(3),
    dwdot_i    = RealVectorX::Zero(3),
    tmp        = RealVectorX::Zero(6);
    
    ComplexVectorX
    vec2_n1    = ComplexVectorX::Zero(n1),
    vec2_n2    = ComplexVectorX::Zero(n2),
    flux       = ComplexVectorX::Zero(n1),
    tmp_c      = ComplexVectorX::Zero(6),
    Dw_i       = ComplexVectorX::Zero(3),
    Dni        = ComplexVectorX::Zero(3);
    
    RealMatrixX
    mat1_n1n1  = RealMatrixX::Zero( n1, n1),
    mat2_n1n2  = RealMatrixX::Zero( n1, n2),
    mat3_n2n2  = RealMatrixX::Zero( n2, n2);
    
    libMesh::Point pt;
    MAST::FEMOperatorMatrix Bmat;
    
    // create objects to calculate the primitive solution, flux, and Jacobian
    MAST::PrimitiveSolution                            primitive_sol;
    MAST::SmallPerturbationPrimitiveSolution<Complex>  sd_primitive_sol;
    
    Complex
    iota            (0., 1.),
    Dvi_ni_freq_indep   = 0.,
    Dvi_ni_freq_dep     = 0.;
    
    Real
    omega               = 0.,
    domega              = 0.,
    b_V                 = 0.;
    
    
    // get the surface motion object from the boundary condition object
    MAST::FieldFunction<RealVectorX>
    *displ   = nullptr;
    MAST::NormalRotationFunctionBase<RealVectorX>
    *n_rot = nullptr;
    
    MAST::FieldFunction<ComplexVectorX>
    *displ_perturb   = nullptr;
    MAST::NormalRotationFunctionBase<ComplexVectorX>
    *n_rot_perturb = nullptr;
    
    
    if (bc.contains("displacement")) {
        
        displ = &bc.get<MAST::FieldFunction<RealVectorX> >("displacement");

        libmesh_assert( bc.contains("normal_rotation"));
        
        MAST::FieldFunction<RealVectorX>&
        tmp = bc.get<MAST::FieldFunction<RealVectorX> >("normal_rotation");
        n_rot = dynamic_cast<MAST::NormalRotationFunctionBase<RealVectorX>*>(&tmp);
    }
    

    if (bc.contains("frequency_domain_displacement")) {
        
        displ_perturb = &bc.get<MAST::FieldFunction<ComplexVectorX> >("frequency_domain_displacement");
        
        libmesh_assert( bc.contains("frequency_domain_normal_rotation"));
        
        MAST::FieldFunction<ComplexVectorX>&
        tmp = bc.get<MAST::FieldFunction<ComplexVectorX> >("frequency_domain_normal_rotation");
        n_rot_perturb = dynamic_cast<MAST::NormalRotationFunctionBase<ComplexVectorX>*>(&tmp);
    }

    
    (*freq)(omega);
    freq->derivative(p, domega);
    freq->nondimensionalizing_factor(b_V);
    
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        // initialize the Bmat operator for this term
        _initialize_fem_interpolation_operator(qp, dim, *fe, Bmat);
        Bmat.right_multiply(vec1_n1,         _sol);  // conservative sol
        Bmat.right_multiply(vec2_n1, _complex_sol);  // perturbation sol
        
        // initialize the primitive solution
        primitive_sol.zero();
        primitive_sol.init(dim,
                           vec1_n1,
                           flight_condition->gas_property.cp,
                           flight_condition->gas_property.cv,
                           if_viscous());
        
        // initialize the small-disturbance primitive sol
        sd_primitive_sol.zero();
        sd_primitive_sol.init(primitive_sol, vec2_n1, if_viscous());
        
        ////////////////////////////////////////////////////////////
        //   Calculation of the surface velocity term. For a
        //   steady-state system,
        //        vi (ni + dni)  =  wdot_i (ni + dni)
        //  or    vi ni = wdot_i (ni + dni) - vi dni
        //
        //    The first order perturbation, D(.), of this is
        //        (vi + Dvi) (ni + dni + Dni) =
        //          (wdot_i + Dwdot_i) (ni + dni + Dni)
        //  or    Dvi ni = Dwdot_i (ni + dni + Dni) +
        //                  wdot_i (ni + dni + Dni) -
        //                  vi (ni + dni + Dni) -
        //                  Dvi (dni + Dni)
        //     Neglecting HOT, this becomes
        //        Dvi ni = Dwdot_i (ni + dni) +
        //                  wdot_i (ni + dni + Dni) -
        //                  vi (ni + dni + Dni) -
        //                  Dvi dni
        //
        ////////////////////////////////////////////////////////////
        
        // copy the surface normal
        for (unsigned int i_dim=0; i_dim<dim; i_dim++)
            ni(i_dim) = normals[qp](i_dim);
        
        // now check if the surface deformation is defined and
        // needs to be applied through transpiration boundary
        // condition
        primitive_sol.get_uvec(uvec);
        
        if (displ) {
            //////////////////////////////////////////////////////////////
            // contribution from the base-flow boundary condition
            //////////////////////////////////////////////////////////////
            (*displ)  (qpoint[qp], _time, tmp);
            dwdot_i = tmp.topRows(3);
            (*n_rot)(qpoint[qp], normals[qp], _time, dni);
        }
        
        if (displ_perturb) {
            
            //////////////////////////////////////////////////////////////
            // contribution from the small-disturbance boundary condition
            //////////////////////////////////////////////////////////////
            (*displ_perturb)(qpoint[qp], _time, tmp_c);
            Dw_i = tmp_c.topRows(3);
            (*n_rot_perturb)(qpoint[qp], normals[qp], _time, Dni);
        }
        
        
        Dvi_ni_freq_dep   = Dw_i.dot(ni+dni) * iota * domega;
        
        flux.setZero();
        flux             += Dvi_ni_freq_dep * vec1_n1;              // Dvi_ni cons_flux
        flux(n1-1)       += Dvi_ni_freq_dep * primitive_sol.p;      // Dvi_ni {0,0,0,0,p}
        
        
        Bmat.vector_mult_transpose(vec2_n2, flux);
        f += JxW[qp] * vec2_n2;
        
        if ( request_jacobian ) {
            libmesh_assert(false); // jac sens not implemented
            
        }
    }
    
    return request_jacobian;
}




