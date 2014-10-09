
// MAST includes
#include "heat_conduction/heat_conduction_elem_base.h"
#include "numerics/fem_operator_matrix.h"
#include "base/system_initialization.h"
#include "base/field_function_base.h"
#include "base/parameter.h"
#include "base/boundary_condition_base.h"
#include "property_cards/element_property_card_base.h"
#include "property_cards/element_property_card_1D.h"
#include "mesh/local_elem_base.h"
#include "mesh/local_1d_elem.h"
#include "mesh/local_2d_elem.h"
#include "mesh/local_3d_elem.h"


MAST::HeatConductionElementBase::HeatConductionElementBase(MAST::SystemInitialization& sys,
                                                           const libMesh::Elem& elem,
                                                           const MAST::ElementPropertyCardBase& p):
MAST::ElementBase(sys, elem),
_property(p) {

    MAST::LocalElemBase* rval = NULL;
    
    switch (elem.dim()) {
        case 1: {
            const MAST::ElementPropertyCard1D& p_1d =
            dynamic_cast<const MAST::ElementPropertyCard1D&>(p);
            rval = new MAST::Local1DElem(elem, p_1d.y_vector());
        }
            break;
            
        case 2:
            rval = new MAST::Local2DElem(elem);
            break;
            
        case 3:
            rval = new MAST::Local3DElem(elem);
            break;
            
        default:
            // should not get here.
            libmesh_error();
            break;
    }
    
    _local_elem.reset(rval);
    
    // now initialize the finite element data structures
    _init_fe_and_qrule(_local_elem->local_elem());
}



MAST::HeatConductionElementBase::~HeatConductionElementBase() {
    
}




bool
MAST::HeatConductionElementBase::internal_residual (bool request_jacobian,
                                                    RealVectorX& f,
                                                    RealMatrixX& jac) {
    const std::vector<Real>& JxW           = _fe->get_JxW();
    const std::vector<libMesh::Point>& xyz = _fe->get_xyz();
    const unsigned int
    n_phi  = _fe->n_shape_functions(),
    dim    = _elem.dim();
    
    RealMatrixX
    material_mat(dim, dim),
    mat_n2n2(n_phi, n_phi);
    RealVectorX
    vec1(1),
    vec2_n2(n_phi),
    flux(dim);
    
    std::auto_ptr<MAST::FieldFunction<RealMatrixX> > conductance =
    _property.thermal_conductance_matrix(*this);
    
    libMesh::Point p;
    std::vector<MAST::FEMOperatorMatrix> dBmat(dim);

    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        _local_elem->global_coordinates_location(xyz[qp], p);
        
        (*conductance)(p, _time, material_mat);
        
        _initialize_flux_fem_operator(qp, dBmat);

        // calculate the flux for each dimension and add its weighted
        // component to the residual
        flux.setZero();
        for (unsigned int j=0; j<dim; j++) {
            dBmat[j].right_multiply(vec1, _sol);
            
            for (unsigned int i=0; i<dim; i++)
                flux(i) += vec1(0) * material_mat(i,j);
        }

        // now add to the residual vector
        for (unsigned int i=0; i<dim; i++) {
            vec1(0)  = flux(i);
            dBmat[i].vector_mult_transpose(vec2_n2, vec1);
            f += JxW[qp] * vec2_n2;
        }

        
        if (request_jacobian) {
            
            // TODO: add Jacobian contribution from temperature dependent property
            for (unsigned int i=0; i<dim; i++)
                for (unsigned int j=0; j<dim; j++) {
                    
                    dBmat[i].right_multiply_transpose(mat_n2n2, dBmat[j]);
                    jac += JxW[qp] * material_mat(i,j) * mat_n2n2;
                }
        }
    }
    
    return request_jacobian;
}





bool
MAST::HeatConductionElementBase::velocity_residual (bool request_jacobian,
                                                    RealVectorX& f,
                                                    RealMatrixX& jac) {
    MAST::FEMOperatorMatrix Bmat;
    
    const std::vector<Real>& JxW                 = _fe->get_JxW();
    const std::vector<libMesh::Point>& xyz       = _fe->get_xyz();
    
    const unsigned int
    n_phi      = _fe->n_shape_functions(),
    dim        = _elem.dim();
    
    RealMatrixX
    material_mat(dim, dim),
    mat_n2n2(n_phi, n_phi);
    RealVectorX
    vec1(1),
    vec2_n2(n_phi);
    
    std::auto_ptr<MAST::FieldFunction<RealMatrixX> > capacitance =
    _property.thermal_capacitance_matrix(*this);
    
    libMesh::Point p;
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        _local_elem->global_coordinates_location(xyz[qp], p);
        
        (*capacitance)(p, _time, material_mat);
        
        _initialize_mass_fem_operator(qp, Bmat);

        Bmat.right_multiply(vec1, _vel);
        Bmat.vector_mult_transpose(vec2_n2, vec1);
        
        f      += JxW[qp] * material_mat(0,0) * vec2_n2;
        
        if (request_jacobian) {
            
            // TODO: add Jacobian contribution from temperature dependent property
            
            Bmat.right_multiply_transpose(mat_n2n2, Bmat);
            jac += JxW[qp] * material_mat(0,0) * mat_n2n2;
        }
        
    }
    
    return request_jacobian;
}




bool
MAST::HeatConductionElementBase::
side_external_residual (bool request_jacobian,
                        RealVectorX& f,
                        RealMatrixX& jac,
                        std::multimap<libMesh::boundary_id_type, MAST::BoundaryConditionBase*>& bc) {
    
    typedef std::multimap<libMesh::boundary_id_type, MAST::BoundaryConditionBase*> maptype;
    
    // iterate over the boundary ids given in the provided force map
    std::pair<maptype::const_iterator, maptype::const_iterator> it;
    
    const libMesh::BoundaryInfo& binfo = *_system.system().get_mesh().boundary_info;
    
    // for each boundary id, check if any of the sides on the element
    // has the associated boundary
    bool calculate_jac = false;
    
    for (unsigned short int n=0; n<_elem.n_sides(); n++) {
        
        // if no boundary ids have been specified for the side, then
        // move to the next side.
        if (!binfo.n_boundary_ids(&_elem, n))
            continue;
        
        // prepare the side finite element
        std::auto_ptr<libMesh::FEBase> fe;
        std::auto_ptr<libMesh::QBase> qrule;
        _get_side_fe_and_qrule(_get_elem_for_quadrature(), n, fe, qrule);
        
        // check to see if any of the specified boundary ids has a boundary
        // condition associated with them
        std::vector<libMesh::boundary_id_type> bc_ids = binfo.boundary_ids(&_elem, n);
        std::vector<libMesh::boundary_id_type>::const_iterator bc_it = bc_ids.begin();
        
        for ( ; bc_it != bc_ids.end(); bc_it++) {
            
            // find the loads on this boundary and evaluate the f and jac
            it = bc.equal_range(*bc_it);
            
            for ( ; it.first != it.second; it.first++) {
                
                // apply all the types of loading
                switch (it.first->second->type()) {
                    case MAST::HEAT_FLUX:
                        calculate_jac = (calculate_jac ||
                                         surface_flux_residual(request_jacobian,
                                                               f, jac,
                                                               *it.first->second,
                                                               *fe));
                        break;
                        
                    case MAST::CONVECTION_HEAT_FLUX:
                        calculate_jac = (calculate_jac ||
                                         surface_convection_residual(request_jacobian,
                                                                     f, jac,
                                                                     *it.first->second,
                                                                     *fe));
                        break;
                        
                    case MAST::SURFACE_RADIATION_HEAT_FLUX:
                        calculate_jac = (calculate_jac ||
                                         surface_radiation_residual(request_jacobian,
                                                                    f, jac,
                                                                    *it.first->second,
                                                                    *fe));
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
    }
    return (request_jacobian && calculate_jac);
}





bool
MAST::HeatConductionElementBase::
volume_external_residual (bool request_jacobian,
                          RealVectorX& f,
                          RealMatrixX& jac,
                          std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*>& bc) {
    
    typedef std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*> maptype;
    
    // iterate over the boundary ids given in the provided force map
    std::pair<maptype::const_iterator, maptype::const_iterator> it;
    
    // for each boundary id, check if any of the sides on the element
    // has the associated boundary
    bool calculate_jac = false;
    
    libMesh::subdomain_id_type sid = _elem.subdomain_id();
    // find the loads on this boundary and evaluate the f and jac
    it =bc.equal_range(sid);
    
    for ( ; it.first != it.second; it.first++) {
        // apply all the types of loading
        switch (it.first->second->type()) {
                
            case MAST::HEAT_FLUX:
                calculate_jac = (calculate_jac ||
                                 surface_flux_residual(request_jacobian,
                                                       f, jac,
                                                       *it.first->second,
                                                       *_fe));
                break;
            case MAST::CONVECTION_HEAT_FLUX:
                calculate_jac = (calculate_jac ||
                                 surface_convection_residual(request_jacobian,
                                                             f, jac,
                                                             *it.first->second,
                                                             *_fe));
                break;
                
            case MAST::SURFACE_RADIATION_HEAT_FLUX:
                calculate_jac = (calculate_jac ||
                                 surface_radiation_residual(request_jacobian,
                                                            f, jac,
                                                            *it.first->second,
                                                            *_fe));
                break;
                
            case MAST::HEAT_SOURCE:
                calculate_jac = (calculate_jac ||
                                 volume_heat_source_residual(request_jacobian,
                                                             f, jac,
                                                             *it.first->second,
                                                             *_fe));
                break;
                
            default:
                // not implemented yet
                libmesh_error();
                break;
        }
    }
    
    return (request_jacobian && calculate_jac);
}



bool
MAST::HeatConductionElementBase::
side_external_residual_sensitivity (bool request_jacobian,
                                    RealVectorX& f,
                                    RealMatrixX& jac,
                                    std::multimap<libMesh::boundary_id_type, MAST::BoundaryConditionBase*>& bc) {
    
    
    typedef std::multimap<libMesh::boundary_id_type, MAST::BoundaryConditionBase*> maptype;
    
    // iterate over the boundary ids given in the provided force map
    std::pair<maptype::const_iterator, maptype::const_iterator> it;
    
    const libMesh::BoundaryInfo& binfo = *_system.system().get_mesh().boundary_info;
    
    // for each boundary id, check if any of the sides on the element
    // has the associated boundary
    bool calculate_jac = false;
    
    for (unsigned short int n=0; n<_elem.n_sides(); n++) {
        
        // if no boundary ids have been specified for the side, then
        // move to the next side.
        if (!binfo.n_boundary_ids(&_elem, n))
            continue;
        
        // prepare the side finite element
        std::auto_ptr<libMesh::FEBase> fe;
        std::auto_ptr<libMesh::QBase> qrule;
        _get_side_fe_and_qrule(_get_elem_for_quadrature(), n, fe, qrule);
        
        // check to see if any of the specified boundary ids has a boundary
        // condition associated with them
        std::vector<libMesh::boundary_id_type> bc_ids = binfo.boundary_ids(&_elem, n);
        std::vector<libMesh::boundary_id_type>::const_iterator bc_it = bc_ids.begin();
        
        for ( ; bc_it != bc_ids.end(); bc_it++) {
            
            // find the loads on this boundary and evaluate the f and jac
            it = bc.equal_range(*bc_it);
            
            for ( ; it.first != it.second; it.first++) {
                
                // apply all the types of loading
                switch (it.first->second->type()) {
                    case MAST::HEAT_FLUX:
                        calculate_jac = (calculate_jac ||
                                         surface_flux_residual_sensitivity(request_jacobian,
                                                                           f, jac,
                                                                           *it.first->second,
                                                                           *fe));
                        break;
                        
                    case MAST::CONVECTION_HEAT_FLUX:
                        calculate_jac = (calculate_jac ||
                                         surface_convection_residual_sensitivity(request_jacobian,
                                                                                 f, jac,
                                                                                 *it.first->second,
                                                                                 *fe));
                        break;
                        
                    case MAST::SURFACE_RADIATION_HEAT_FLUX:
                        calculate_jac = (calculate_jac ||
                                         surface_radiation_residual_sensitivity(request_jacobian,
                                                                                f, jac,
                                                                                *it.first->second,
                                                                                *fe));
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
    }
    return (request_jacobian && calculate_jac);
}




bool
MAST::HeatConductionElementBase::
volume_external_residual_sensitivity (bool request_jacobian,
                                      RealVectorX& f,
                                      RealMatrixX& jac,
                                      std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*>& bc) {
    
    typedef std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*> maptype;
    
    // iterate over the boundary ids given in the provided force map
    std::pair<maptype::const_iterator, maptype::const_iterator> it;
    
    // for each boundary id, check if any of the sides on the element
    // has the associated boundary
    bool calculate_jac = false;
    
    libMesh::subdomain_id_type sid = _elem.subdomain_id();
    // find the loads on this boundary and evaluate the f and jac
    it =bc.equal_range(sid);
    
    for ( ; it.first != it.second; it.first++) {
        // apply all the types of loading
        switch (it.first->second->type()) {
                
            case MAST::HEAT_FLUX:
                calculate_jac = (calculate_jac ||
                                 surface_flux_residual_sensitivity(request_jacobian,
                                                                   f, jac,
                                                                   *it.first->second,
                                                                   *_fe));
                break;
            case MAST::CONVECTION_HEAT_FLUX:
                calculate_jac = (calculate_jac ||
                                 surface_convection_residual_sensitivity(request_jacobian,
                                                                         f, jac,
                                                                         *it.first->second,
                                                                         *_fe));
                break;
                
            case MAST::SURFACE_RADIATION_HEAT_FLUX:
                calculate_jac = (calculate_jac ||
                                 surface_radiation_residual_sensitivity(request_jacobian,
                                                                        f, jac,
                                                                        *it.first->second,
                                                                        *_fe));
                break;
                
            case MAST::HEAT_SOURCE:
                calculate_jac = (calculate_jac ||
                                 volume_heat_source_residual_sensitivity(request_jacobian,
                                                                         f, jac,
                                                                         *it.first->second,
                                                                         *_fe));
                break;
                
            default:
                // not implemented yet
                libmesh_error();
                break;
        }
    }
    
    return (request_jacobian && calculate_jac);
}




bool
MAST::HeatConductionElementBase::
internal_residual_sensitivity (bool request_jacobian,
                               RealVectorX& f,
                               RealMatrixX& jac) {
    
    return request_jacobian;
}



bool
MAST::HeatConductionElementBase::
velocity_residual_sensitivity (bool request_jacobian,
                               RealVectorX& f,
                               RealMatrixX& jac) {
    
    return request_jacobian;
}








bool
MAST::HeatConductionElementBase::
surface_flux_residual(bool request_jacobian,
                      RealVectorX& f,
                      RealMatrixX& jac,
                      MAST::BoundaryConditionBase& p,
                      const libMesh::FEBase& fe) {
    
    libmesh_assert(_elem.dim() < 3); // only applicable for lower dimensional elements
    
    
    // get the function from this boundary condition
    const MAST::FieldFunction<Real>& func =
    p.get<MAST::FieldFunction<Real> >("heat_flux");
    
    
    const std::vector<Real> &JxW                 = fe.get_JxW();
    const std::vector<libMesh::Point>& qpoint    = fe.get_xyz();
    const std::vector<std::vector<Real> >& phi   = fe.get_phi();
    const unsigned int n_phi                     = (unsigned int)phi.size();
    
    RealVectorX phi_vec   (n_phi);
    libMesh::Point pt;
    Real  flux;
    
    for (unsigned int qp=0; qp<qpoint.size(); qp++) {
        
        _local_elem->global_coordinates_location                (qpoint[qp], pt);
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        // get the value of flux = q_i . n_i
        func(pt, _time, flux);
        
        f   +=  JxW[qp] * phi_vec * flux;
    }
    
    // calculation of the load vector is independent of solution
    return false;
}





bool
MAST::HeatConductionElementBase::
surface_flux_residual_sensitivity(bool request_jacobian,
                                  RealVectorX& f,
                                  RealMatrixX& jac,
                                  MAST::BoundaryConditionBase& p,
                                  const libMesh::FEBase& fe) {
    
    return false;
}




bool
MAST::HeatConductionElementBase::
surface_convection_residual(bool request_jacobian,
                            RealVectorX& f,
                            RealMatrixX& jac,
                            MAST::BoundaryConditionBase& p,
                            const libMesh::FEBase& fe) {
    
    // get the function from this boundary condition
    const MAST::FieldFunction<Real>
    &coeff = p.get<MAST::FieldFunction<Real> >("convection_coeff"),
    &T_amb = p.get<MAST::FieldFunction<Real> >("ambient_temperature");
    
    const std::vector<Real> &JxW               = fe.get_JxW();
    const std::vector<libMesh::Point>& qpoint  = fe.get_xyz();
    const std::vector<std::vector<Real> >& phi = fe.get_phi();
    const unsigned int n_phi                   = (unsigned int)phi.size();
    
    
    RealVectorX  phi_vec  (n_phi);
    RealMatrixX  mat      (n_phi, n_phi);
    Real temp, amb_temp, h_coeff;
    libMesh::Point pt;
    MAST::FEMOperatorMatrix Bmat;
    
    
    for (unsigned int qp=0; qp<qpoint.size(); qp++) {
        
        _local_elem->global_coordinates_location (qpoint[qp], pt);
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        // value of flux
        coeff(pt, _time, h_coeff);
        T_amb(pt, _time, amb_temp);
        temp  = phi_vec.dot(_sol);
        
        // normal flux is given as:
        // qi_ni = h_coeff * (T - T_amb)
        //
        f   += JxW[qp] * phi_vec * h_coeff * (temp - amb_temp);
        
        if (request_jacobian) {
            
            Bmat.reinit(1, phi_vec);
            Bmat.right_multiply_transpose(mat, Bmat);
            jac += JxW[qp] * mat * h_coeff;
        }
    }
    
    return request_jacobian;
}





bool
MAST::HeatConductionElementBase::
surface_convection_residual_sensitivity(bool request_jacobian,
                                        RealVectorX& f,
                                        RealMatrixX& jac,
                                        MAST::BoundaryConditionBase& p,
                                        const libMesh::FEBase& fe) {
    
    return false;
}





bool
MAST::HeatConductionElementBase::
surface_radiation_residual(bool request_jacobian,
                           RealVectorX& f,
                           RealMatrixX& jac,
                           MAST::BoundaryConditionBase& p,
                           const libMesh::FEBase& fe) {
    
    // get the function from this boundary condition
    const MAST::FieldFunction<Real>
    &emissivity = p.get<MAST::FieldFunction<Real> >("emissivity"),
    &T_amb      = p.get<MAST::FieldFunction<Real> >("ambient_temperature");
    const MAST::Parameter
    &sb_const   = p.get<MAST::Parameter>("stefan_bolzmann_constant");
    
    
    const std::vector<Real> &JxW               = fe.get_JxW();
    const std::vector<libMesh::Point>& qpoint  = fe.get_xyz();
    const std::vector<std::vector<Real> >& phi = fe.get_phi();
    const unsigned int n_phi                   = (unsigned int)phi.size();
    
    RealVectorX phi_vec  (n_phi);
    RealMatrixX mat      (n_phi, n_phi);
    const Real sbc = sb_const();
    Real temp, amb_temp, emiss;
    libMesh::Point pt;
    MAST::FEMOperatorMatrix Bmat;
    
    for (unsigned int qp=0; qp<qpoint.size(); qp++) {
        
        _local_elem->global_coordinates_location (qpoint[qp], pt);
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        // value of flux
        emissivity(pt, _time, emiss);
        T_amb     (pt, _time, amb_temp);
        temp  = phi_vec.dot(_sol);
        
        f   += JxW[qp] * phi_vec * sbc * emiss * pow(temp-amb_temp, 4.);
        
        if (request_jacobian) {
            
            Bmat.reinit(1, phi_vec);
            Bmat.right_multiply_transpose(mat, Bmat);
            jac +=  JxW[qp] * mat * sbc * emiss * 4. * pow(temp-amb_temp, 3.);
        }
    }
    
    
    return request_jacobian;
}






bool
MAST::HeatConductionElementBase::
surface_radiation_residual_sensitivity(bool request_jacobian,
                                       RealVectorX& f,
                                       RealMatrixX& jac,
                                       MAST::BoundaryConditionBase& p,
                                       const libMesh::FEBase& fe) {
    
    return true;
}





bool
MAST::HeatConductionElementBase::
volume_heat_source_residual(bool request_jacobian,
                            RealVectorX& f,
                            RealMatrixX& jac,
                            MAST::BoundaryConditionBase& p,
                            const libMesh::FEBase& fe) {
    
    return false;
}



bool
MAST::HeatConductionElementBase::
volume_heat_source_residual_sensitivity(bool request_jacobian,
                                        RealVectorX& f,
                                        RealMatrixX& jac,
                                        MAST::BoundaryConditionBase& p,
                                        const libMesh::FEBase& fe) {
    
    return false;
}



void
MAST::HeatConductionElementBase::
_initialize_mass_fem_operator(const unsigned int qp,
                              MAST::FEMOperatorMatrix& Bmat) {
    
    const std::vector<std::vector<Real> >& phi_fe = _fe->get_phi();
    
    const unsigned int n_phi = (unsigned int)phi_fe.size();
    
    RealVectorX phi(n_phi);
    
    // shape function values
    // N
    for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
        phi(i_nd) = phi_fe[i_nd][qp];
    
    Bmat.reinit(1, phi);
}




void
MAST::HeatConductionElementBase::
_initialize_flux_fem_operator(const unsigned int qp,
                              std::vector<MAST::FEMOperatorMatrix>& dBmat) {
    
    const std::vector<std::vector<libMesh::RealVectorValue> >& dphi = _fe->get_dphi();
    
    const unsigned int n_phi = (unsigned int)dphi.size();
    RealVectorX phi(n_phi);
    
    // now set the shape function values
    // dN/dx
    for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
        phi(i_nd) = dphi[i_nd][qp](0);
    
    dBmat[0].reinit(1, phi); // dT/dx
    
    // only for 2D and 3D elements
    if (_elem.dim() > 1) {
        
        // dN/dy
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi(i_nd) = dphi[i_nd][qp](1);
        dBmat[1].reinit(1, phi); //  dT/dy
        
        // only for 3D elements
        if (_elem.dim() > 2) {
            
            // dN/dz
            for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
                phi(i_nd) = dphi[i_nd][qp](2);
            dBmat[2].reinit(1, phi); //  dT/dz
        }
    }
}


