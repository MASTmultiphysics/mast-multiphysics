
// MAST includes
#include "heat_conduction/heat_conduction_elem_base.h"
#include "numerics/fem_operator_matrix.h"
#include "base/system_initialization.h"
#include "base/field_function_base.h"
#include "base/boundary_condition_base.h"
#include "property_cards/element_property_card_base.h"
#include "mesh/local_elem_base.h"


MAST::HeatConductionElementBase::HeatConductionElementBase(MAST::SystemInitialization& sys,
                                                           const libMesh::Elem& elem,
                                                           MAST::ElementPropertyCardBase& p):
MAST::ElementBase(sys, elem),
_property(p) {
    
}



MAST::HeatConductionElementBase::~HeatConductionElementBase() {
    
}




bool
MAST::HeatConductionElementBase::internal_residual (bool request_jacobian,
                                                    RealVectorX& f,
                                                    RealMatrixX& jac) {
    MAST::FEMOperatorMatrix Bmat;
    
    const std::vector<Real>& JxW           = _fe->get_JxW();
    const std::vector<libMesh::Point>& xyz = _fe->get_xyz();
    const unsigned int
    n_phi  = _fe->n_shape_functions(),
    dim    = _elem.dim();
    
    RealMatrixX
    material_mat(dim, dim),
    mat1_n1n2(dim, n_phi),
    mat2_n2n2(dim, n_phi);
    RealVectorX
    phi_vec(n_phi),
    vec1_n1(dim),
    vec2_n2(n_phi);
    
    std::auto_ptr<MAST::FieldFunction<RealMatrixX> > conductance =
    _property.thermal_conductance_matrix(*this);
    
    libMesh::Point p;
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        _local_elem->global_coordinates_location(xyz[qp], p);
        
        (*conductance)(p, _time, material_mat);
        
        _initialize_fem_operator(qp, false, Bmat);
        
        Bmat.left_multiply(mat1_n1n2, material_mat);
        
        vec1_n1 = mat1_n1n2 * _sol;
        Bmat.vector_mult_transpose(vec2_n2, vec1_n1);
        f += JxW[qp] * vec2_n2;
        
        if (request_jacobian) {
            
            // TODO: add Jacobian contribution from temperature dependent property
            
            Bmat.right_multiply_transpose(mat2_n2n2, mat1_n1n2);
            jac += JxW[qp] * mat2_n2n2;
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
    mat1_n1n2(dim, n_phi),
    mat2_n2n2(dim, n_phi);
    RealVectorX
    phi_vec(n_phi),
    vec1_n1(dim),
    vec2_n2(n_phi);
    
    std::auto_ptr<MAST::FieldFunction<RealMatrixX> > capacitance =
    _property.thermal_capacitance_matrix(*this);
    
    libMesh::Point p;
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        _local_elem->global_coordinates_location(xyz[qp], p);
        
        (*capacitance)(p, _time, material_mat);
        
        _initialize_fem_operator(qp, true, Bmat);
        
        Bmat.left_multiply(mat1_n1n2, material_mat);
        
        vec1_n1 = mat1_n1n2 * _vel;
        Bmat.vector_mult_transpose(vec2_n2, vec1_n1);
        f      += JxW[qp] * vec2_n2;
        
        if (request_jacobian) {
            
            // TODO: add Jacobian contribution from temperature dependent property
            
            Bmat.right_multiply_transpose(mat2_n2n2, mat1_n1n2);
            jac += JxW[qp] * mat2_n2n2;
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
                                                               n,
                                                               *it.first->second));
                        break;
                        
                    case MAST::CONVECTION_HEAT_FLUX:
                        calculate_jac = (calculate_jac ||
                                         surface_convection_residual(request_jacobian,
                                                                     f, jac,
                                                                     n,
                                                                     *it.first->second));
                        break;
                        
                    case MAST::SURFACE_RADIATION_HEAT_FLUX:
                        calculate_jac = (calculate_jac ||
                                         surface_radiation_residual(request_jacobian,
                                                                    f, jac,
                                                                    n,
                                                                    *it.first->second));
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
                                                       *it.first->second));
                break;
            case MAST::CONVECTION_HEAT_FLUX:
                calculate_jac = (calculate_jac ||
                                 surface_convection_residual(request_jacobian,
                                                             f, jac,
                                                             *it.first->second));
                break;
                
            case MAST::SURFACE_RADIATION_HEAT_FLUX:
                calculate_jac = (calculate_jac ||
                                 surface_radiation_residual(request_jacobian,
                                                            f, jac,
                                                            *it.first->second));
                break;

            case MAST::HEAT_SOURCE:
                calculate_jac = (calculate_jac ||
                                 volume_heat_source_residual(request_jacobian,
                                                             f, jac,
                                                             *it.first->second));
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
                                                                           n,
                                                                           *it.first->second));
                        break;
                        
                    case MAST::CONVECTION_HEAT_FLUX:
                        calculate_jac = (calculate_jac ||
                                         surface_convection_residual_sensitivity(request_jacobian,
                                                                                 f, jac,
                                                                                 n,
                                                                                 *it.first->second));
                        break;
                        
                    case MAST::SURFACE_RADIATION_HEAT_FLUX:
                        calculate_jac = (calculate_jac ||
                                         surface_radiation_residual_sensitivity(request_jacobian,
                                                                                f, jac,
                                                                                n,
                                                                                *it.first->second));
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
                                                       *it.first->second));
                break;
            case MAST::CONVECTION_HEAT_FLUX:
                calculate_jac = (calculate_jac ||
                                 surface_convection_residual_sensitivity(request_jacobian,
                                                                         f, jac,
                                                                         *it.first->second));
                break;
                
            case MAST::SURFACE_RADIATION_HEAT_FLUX:
                calculate_jac = (calculate_jac ||
                                 surface_radiation_residual_sensitivity(request_jacobian,
                                                                        f, jac,
                                                                        *it.first->second));
                break;
                
            case MAST::HEAT_SOURCE:
                calculate_jac = (calculate_jac ||
                                 volume_heat_source_residual_sensitivity(request_jacobian,
                                                                         f, jac,
                                                                         *it.first->second));
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
                               RealMatrixX& jac,
                               bool if_ignore_ho_jac) {
    
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
                      const unsigned int side,
                      MAST::BoundaryConditionBase& p) {
    
    // get the function from this boundary condition
    const MAST::FieldFunction<RealVector3>& func =
    p.get<MAST::FieldFunction<RealVector3> >("heat_flux");
    
    std::auto_ptr<libMesh::FEBase> fe;
    std::auto_ptr<libMesh::QBase> qrule;
    _get_side_fe_and_qrule(_get_elem_for_quadrature(), side, fe, qrule);
    
    
    const std::vector<Real> &JxW               = fe->get_JxW();
    const std::vector<libMesh::Point>& qpoint  = fe->get_xyz();
    const std::vector<std::vector<Real> >& phi = fe->get_phi();
    const unsigned int
    n_phi      = (unsigned int)phi.size();
    
    // boundary normals
    const std::vector<libMesh::Point>& face_normals = fe->get_normals();
    
    RealVectorX
    phi_vec  (n_phi);
    RealVector3
    flux     (3),
    normal   (3);
    libMesh::Point pt;
    
    for (unsigned int qp=0; qp<qpoint.size(); qp++) {
        
        _local_elem->global_coordinates_location(qpoint[qp], pt);
        _local_elem->global_coordinates_normal(face_normals[qp], normal);
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        // value of flux
        func(pt, _time, flux);
        
        f   += JxW[qp] * phi_vec * flux.dot(normal);
    }
    
    // calculation of the load vector is independent of solution
    return false;
}




bool
MAST::HeatConductionElementBase::
surface_flux_residual(bool request_jacobian,
                      RealVectorX& f,
                      RealMatrixX& jac,
                      MAST::BoundaryConditionBase& p) {
    
    libmesh_assert(_elem.dim() < 3); // only applicable for lower dimensional elements
    
    
    // get the function from this boundary condition
    const MAST::FieldFunction<RealVector3>& func =
    p.get<MAST::FieldFunction<RealVector3> >("heat_flux");
    
    
    const std::vector<Real> &JxW                 = _fe->get_JxW();
    const std::vector<libMesh::Point>& qpoint    = _fe->get_xyz();
    const std::vector<std::vector<Real> >& phi   = _fe->get_phi();
    const unsigned int
    n_phi        = (unsigned int)phi.size();
    
    RealVectorX
    phi_vec   (n_phi);
    RealVector3
    flux      (3),
    normal    (3);
    libMesh::Point pt;
    
    for (unsigned int qp=0; qp<qpoint.size(); qp++) {
        
        _local_elem->global_coordinates_location(qpoint[qp], pt);
        _local_elem->global_coordinates_normal(pt, normal);
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        // get pressure value
        func(pt, _time, flux);
        
        f   +=  JxW[qp] * phi_vec * flux.dot(normal);
    }
    
    // calculation of the load vector is independent of solution
    return false;
}




bool
MAST::HeatConductionElementBase::
surface_flux_residual_sensitivity(bool request_jacobian,
                                  RealVectorX& f,
                                  RealMatrixX& jac,
                                  const unsigned int side,
                                  MAST::BoundaryConditionBase& p) {
    
    return false;
}




bool
MAST::HeatConductionElementBase::
surface_flux_residual_sensitivity(bool request_jacobian,
                                  RealVectorX& f,
                                  RealMatrixX& jac,
                                  MAST::BoundaryConditionBase& p) {
    
    return false;
}




bool
MAST::HeatConductionElementBase::
surface_convection_residual(bool request_jacobian,
                            RealVectorX& f,
                            RealMatrixX& jac,
                            const unsigned int side,
                            MAST::BoundaryConditionBase& p) {
    
    // get the function from this boundary condition
    const MAST::FieldFunction<RealVector3>& func =
    p.get<MAST::FieldFunction<RealVector3> >("heat_flux");
    
    std::auto_ptr<libMesh::FEBase> fe;
    std::auto_ptr<libMesh::QBase> qrule;
    _get_side_fe_and_qrule(_get_elem_for_quadrature(), side, fe, qrule);
    
    
    const std::vector<Real> &JxW               = fe->get_JxW();
    const std::vector<libMesh::Point>& qpoint  = fe->get_xyz();
    const std::vector<std::vector<Real> >& phi = fe->get_phi();
    const unsigned int
    n_phi      = (unsigned int)phi.size();
    
    // boundary normals
    const std::vector<libMesh::Point>& face_normals = fe->get_normals();
    
    RealVectorX
    phi_vec  (n_phi);
    RealVector3
    flux     (3),
    normal   (3);
    libMesh::Point pt;
    
    for (unsigned int qp=0; qp<qpoint.size(); qp++) {
        
        _local_elem->global_coordinates_location(qpoint[qp], pt);
        _local_elem->global_coordinates_normal(face_normals[qp], normal);
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        // value of flux
        func(pt, _time, flux);
        
        f   += JxW[qp] * phi_vec * flux.dot(normal);
    }
    
    // calculation of the load vector is independent of solution
    return false;
}




bool
MAST::HeatConductionElementBase::
surface_convection_residual(bool request_jacobian,
                            RealVectorX& f,
                            RealMatrixX& jac,
                            MAST::BoundaryConditionBase& p) {
    
    libmesh_assert(_elem.dim() < 3); // only applicable for lower dimensional elements
    
    
    // get the function from this boundary condition
    const MAST::FieldFunction<RealVector3>& func =
    p.get<MAST::FieldFunction<RealVector3> >("heat_flux");
    
    
    const std::vector<Real> &JxW                 = _fe->get_JxW();
    const std::vector<libMesh::Point>& qpoint    = _fe->get_xyz();
    const std::vector<std::vector<Real> >& phi   = _fe->get_phi();
    const unsigned int
    n_phi        = (unsigned int)phi.size();
    
    RealVectorX
    phi_vec   (n_phi);
    RealVector3
    flux      (3),
    normal    (3);
    libMesh::Point pt;
    
    for (unsigned int qp=0; qp<qpoint.size(); qp++) {
        
        _local_elem->global_coordinates_location(qpoint[qp], pt);
        _local_elem->domain_surface_normal_in_global_coordinates(pt, normal);
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        // get pressure value
        func(pt, _time, flux);
        
        f   +=  JxW[qp] * phi_vec * flux.dot(normal);
    }
    
    // calculation of the load vector is independent of solution
    return false;
}




bool
MAST::HeatConductionElementBase::
surface_convection_residual_sensitivity(bool request_jacobian,
                                  RealVectorX& f,
                                  RealMatrixX& jac,
                                  const unsigned int side,
                                  MAST::BoundaryConditionBase& p) {
    
    return false;
}




bool
MAST::HeatConductionElementBase::
surface_convection_residual_sensitivity(bool request_jacobian,
                                        RealVectorX& f,
                                        RealMatrixX& jac,
                                        MAST::BoundaryConditionBase& p) {
    
    return false;
}




bool
MAST::HeatConductionElementBase::
surface_radiation_residual(bool request_jacobian,
                           RealVectorX& f,
                           RealMatrixX& jac,
                           const unsigned int side,
                           MAST::BoundaryConditionBase& p) {
    
    // get the function from this boundary condition
    const MAST::FieldFunction<RealVector3>& func =
    p.get<MAST::FieldFunction<RealVector3> >("heat_flux");
    
    std::auto_ptr<libMesh::FEBase> fe;
    std::auto_ptr<libMesh::QBase> qrule;
    _get_side_fe_and_qrule(_get_elem_for_quadrature(), side, fe, qrule);
    
    
    const std::vector<Real> &JxW               = fe->get_JxW();
    const std::vector<libMesh::Point>& qpoint  = fe->get_xyz();
    const std::vector<std::vector<Real> >& phi = fe->get_phi();
    const unsigned int
    n_phi      = (unsigned int)phi.size();
    
    // boundary normals
    const std::vector<libMesh::Point>& face_normals = fe->get_normals();
    
    RealVectorX
    phi_vec  (n_phi);
    RealVector3
    flux     (3),
    normal   (3);
    libMesh::Point pt;
    
    for (unsigned int qp=0; qp<qpoint.size(); qp++) {
        
        _local_elem->global_coordinates_location(qpoint[qp], pt);
        _local_elem->global_coordinates_normal(face_normals[qp], normal);
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        // value of flux
        func(pt, _time, flux);
        
        f   += JxW[qp] * phi_vec * flux.dot(normal);
    }
    
    // calculation of the load vector is independent of solution
    return true;
}




bool
MAST::HeatConductionElementBase::
surface_radiation_residual(bool request_jacobian,
                           RealVectorX& f,
                           RealMatrixX& jac,
                           MAST::BoundaryConditionBase& p) {
    
    libmesh_assert(_elem.dim() < 3); // only applicable for lower dimensional elements
    
    
    // get the function from this boundary condition
    const MAST::FieldFunction<RealVector3>& func =
    p.get<MAST::FieldFunction<RealVector3> >("heat_flux");
    
    
    const std::vector<Real> &JxW                 = _fe->get_JxW();
    const std::vector<libMesh::Point>& qpoint    = _fe->get_xyz();
    const std::vector<std::vector<Real> >& phi   = _fe->get_phi();
    const unsigned int
    n_phi        = (unsigned int)phi.size();
    
    RealVectorX
    phi_vec   (n_phi);
    RealVector3
    flux      (3),
    normal    (3);
    libMesh::Point pt;
    
    for (unsigned int qp=0; qp<qpoint.size(); qp++) {
        
        _local_elem->global_coordinates_location(qpoint[qp], pt);
        _local_elem->domain_surface_normal_in_global_coordinates(pt, normal);
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        // get pressure value
        func(pt, _time, flux);
        
        f   +=  JxW[qp] * phi_vec * flux.dot(normal);
    }
    
    // calculation of the load vector is independent of solution
    return true;
}




bool
MAST::HeatConductionElementBase::
surface_radiation_residual_sensitivity(bool request_jacobian,
                                       RealVectorX& f,
                                       RealMatrixX& jac,
                                       const unsigned int side,
                                       MAST::BoundaryConditionBase& p) {
    
    return true;
}




bool
MAST::HeatConductionElementBase::
surface_radiation_residual_sensitivity(bool request_jacobian,
                                       RealVectorX& f,
                                       RealMatrixX& jac,
                                       MAST::BoundaryConditionBase& p) {
    
    return true;
}



bool
MAST::HeatConductionElementBase::
volume_heat_source_residual(bool request_jacobian,
                            RealVectorX& f,
                            RealMatrixX& jac,
                            MAST::BoundaryConditionBase& p) {
    
    return false;
}



bool
MAST::HeatConductionElementBase::
volume_heat_source_residual_sensitivity(bool request_jacobian,
                                        RealVectorX& f,
                                        RealMatrixX& jac,
                                        MAST::BoundaryConditionBase& p) {
    
    return false;
}



void
MAST::HeatConductionElementBase::
_initialize_fem_operator(const unsigned int qp,
                         const bool mass,
                         MAST::FEMOperatorMatrix& Bmat) {
    
    if (mass) {
        
        const std::vector<std::vector<Real> >& phi_fe = _fe->get_phi();
        
        const unsigned int n_phi = (unsigned int)phi_fe.size();
        
        RealVectorX phi(n_phi);
        
        
        // shape function values
        // N
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi(i_nd) = phi_fe[i_nd][qp];
        
        for (unsigned int i_dim=0; i_dim<_elem.dim(); i_dim++)
            Bmat.set_shape_function(i_dim, 0, phi); // T
    }
    else {
        
        const std::vector<std::vector<libMesh::RealVectorValue> >& dphi = _fe->get_dphi();
        
        const unsigned int n_phi = (unsigned int)dphi.size();
        RealVectorX phi(n_phi);
        
        
        // now set the shape function values
        // dN/dx
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi(i_nd) = dphi[i_nd][qp](0);
        
        Bmat.set_shape_function(0, 0, phi); // dT/dx
        
        // only for 2D and 3D elements
        if (_elem.dim() > 1) {
            
            // dN/dy
            for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
                phi(i_nd) = dphi[i_nd][qp](1);
            Bmat.set_shape_function(1, 0, phi); //  dT/dy
            
            // only for 3D elements
            if (_elem.dim() > 2) {
                
                // dN/dz
                for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
                    phi(i_nd) = dphi[i_nd][qp](2);
                Bmat.set_shape_function(2, 0, phi); //  dT/dz
            }
        }
    }
}


