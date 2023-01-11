/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2020  Manav Bhatia and MAST authors
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
#include "elasticity/center_of_mass_output.h"
#include "elasticity/structural_element_base.h"
#include "base/assembly_base.h"
#include "base/constant_field_function.h"
#include "base/field_function_base.h"
#include "base/system_initialization.h"
#include "base/nonlinear_system.h"
#include "base/physics_discipline_base.h"
#include "base/boundary_condition_base.h"
#include "base/parameter.h"
#include "property_cards/element_property_card_1D.h"
#include "property_cards/material_property_card_base.h"
#include "level_set/level_set_intersection.h"
#include "level_set/level_set_intersected_elem.h"
#include "mesh/geom_elem.h"
#include "mesh/fe_base.h"

// libMesh includes
#include "libmesh/parallel.h"


MAST::CenterOfMassOutput::CenterOfMassOutput():
MAST::OutputAssemblyElemOperations() {
    _center_of_mass = RealVectorX::Zero(3,1);
    _dcenter_of_mass_dp = RealVectorX::Zero(3,1);
    _first_moment_of_inertia_yz = 0.0;
    _first_moment_of_inertia_xz = 0.0;
    _first_moment_of_inertia_xy = 0.0;
    _dfirst_moment_of_inertia_yz_dp = 0.0;
    _dfirst_moment_of_inertia_xz_dp = 0.0;
    _dfirst_moment_of_inertia_xy_dp = 0.0;
    _mass = 0.0;
    _dmass_dp = 0.0;
}


MAST::CenterOfMassOutput::~CenterOfMassOutput() {

}


void
MAST::CenterOfMassOutput::zero_for_analysis() {
    _center_of_mass.setZero();
    _dcenter_of_mass_dp.setZero();
    _first_moment_of_inertia_yz = 0.0;
    _first_moment_of_inertia_xz = 0.0;
    _first_moment_of_inertia_xy = 0.0;
    _mass = 0.0;

}


void
MAST::CenterOfMassOutput::zero_for_sensitivity() {
    
    _center_of_mass.setZero();
    _dcenter_of_mass_dp.setZero();
}


void
MAST::CenterOfMassOutput::evaluate() {
    
    // make sure that this has not been initialized and calculated for all elems
    libmesh_assert(_physics_elem);

    if (this->if_evaluate_for_element(_physics_elem->elem())) {

        // Get the density field function
        const MAST::ElementPropertyCardBase& section_property = dynamic_cast<const MAST::ElementPropertyCardBase&>(_discipline->get_property_card(_physics_elem->elem()));
        const MAST::FieldFunction<Real>& rho_f  = section_property.get_material().get<MAST::FieldFunction<Real>>("rho");
        // Another method of doing this below
        // MAST::StructuralElementBase& e = dynamic_cast<MAST::StructuralElementBase&>(*_physics_elem);
        // const MAST::FieldFunction<Real>& rho_f = const_cast<MAST::MaterialPropertyCardBase&>(e.elem_property().get_material()).get<MAST::FieldFunction<Real>>("rho");

        // Get the integration points and weights
        // TODO: Currently using the same integration scheme as for element building. Need to allow using user-specified integration scheme
        std::unique_ptr<MAST::FEBase> fe(_physics_elem->elem().init_fe(true, false));
        std::vector<Real> JxW = fe->get_JxW();
        const std::vector<libMesh::Point>& xyz = fe->get_xyz();
        const unsigned int n_qp = (unsigned int)fe->get_qpoints().size();

        // Initialize some variables
        Real rho = 0.0;
        Real Vc = 0.0;
        const MAST::FieldFunction<Real>* Vc_f;

        // Create a unity field function that evaluate to 1 over the entire domain, used for 3D elements
        MAST::Parameter one("one", 1.0);
        MAST::ConstantFieldFunction one_f("one", one);

        // Create a field function to calculate the area of a 1D rectangular element (old MAST)
        class Area1DFieldFunction: public MAST::FieldFunction<Real> {
        public:
            Area1DFieldFunction(const MAST::FieldFunction<Real>& hy, const MAST::FieldFunction<Real>& hz):
            MAST::FieldFunction<Real>("A"), _hy(hy), _hz(hz) {}
            virtual void operator()(const libMesh::Point& p, const Real t, Real& v) const {
                Real hy;
                Real hz;
                _hy(p, t, hy);
                _hz(p, t, hz);
                v = hy * hz;
            }
            virtual void derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t, Real& v) const {
                Real hy, hz, dhy, dhz;
                if (&f == &_hy) {
                    _hz(p, t, hz);
                    v = hz;
                }
                else if (&f == &_hz) {
                    _hy(p, t, hy);
                    v = hy;
                }
                else {
                    _hy(p, t, hy);
                    _hy.derivative(f, p, t, dhy);
                    _hz(p, t, hz);
                    _hz.derivative(f, p, t, dhz);
                    v = dhy * hz + hy * dhz;
                }
            }
        protected:
            const MAST::FieldFunction<Real>& _hy;
            const MAST::FieldFunction<Real>& _hz;
        };
        Area1DFieldFunction A1d_f(section_property.get<MAST::FieldFunction<Real>>("hy"),
                                  section_property.get<MAST::FieldFunction<Real>>("hz"));
        // TODO: Once the 1D cross-section library is implemented in this branch, the class and field function above will no longer be needed.
        //   this is only used for "old" MAST where only rectangular cross sections were supported.

        // To get the total volume for integration, we need to multiply by area for 1D elements and by thickness for 2D elements
        if (section_property.dim() == 1) {
            // Get the area for 1D elements since volume = JxW * A in this case
            // Vc_f = &(section_property.get<MAST::FieldFunction<Real>>("A"));
            Vc_f = &A1d_f;
        }
        else if (section_property.dim() == 2) {
            // Get the thickness for 2D elements since volume = JxW * h in this case
            Vc_f = &(section_property.get<MAST::FieldFunction<Real>>("h"));
        }
        else if (section_property.dim() == 3) {
            // Get unity (1) for 3D element since volume = JxW in this case
            Vc_f = &one_f;
        }

        Real mass_elem = 0.0;
        Real time = 0.0;
        // Iterate over the integration points to perform the numerical integration
        for (unsigned int ind_qp=0; ind_qp<n_qp; ind_qp++) {
            // Get the current time of the system
            time = _system->system().time;

            // Get the density at this integration point
            rho_f(xyz[ind_qp], time, rho);
            // Get the extra volume component at this integration point
            (*Vc_f)(xyz[ind_qp], time, Vc);

            mass_elem = rho * Vc * JxW[ind_qp]; // A value used in multiple calculations below

            _first_moment_of_inertia_yz += mass_elem * xyz[ind_qp](0);
            _first_moment_of_inertia_xz += mass_elem * xyz[ind_qp](1);
            _first_moment_of_inertia_xy += mass_elem * xyz[ind_qp](2);
            _mass += mass_elem;
        }
    }
}


void
MAST::CenterOfMassOutput::evaluate_sensitivity(const MAST::FunctionBase &f) {
    
    // make sure that this has not been initialized ana calculated for all elems
    libmesh_assert(_physics_elem);
    
        if (this->if_evaluate_for_element(_physics_elem->elem())) {

        // Get the density field function
        const MAST::ElementPropertyCardBase& section_property = dynamic_cast<const MAST::ElementPropertyCardBase&>(_discipline->get_property_card(_physics_elem->elem()));
        const MAST::FieldFunction<Real>& rho_f  = section_property.get_material().get<MAST::FieldFunction<Real>>("rho");
        // Another method of doing this below
        // MAST::StructuralElementBase& e = dynamic_cast<MAST::StructuralElementBase&>(*_physics_elem);
        // const MAST::FieldFunction<Real>& rho_f = const_cast<MAST::MaterialPropertyCardBase&>(e.elem_property().get_material()).get<MAST::FieldFunction<Real>>("rho");

        // Get the integration points and weights
        // TODO: Currently using the same integration scheme as for element building. Need to allow using user-specified integration scheme
        std::unique_ptr<MAST::FEBase> fe(_physics_elem->elem().init_fe(true, false));
        std::vector<Real> JxW = fe->get_JxW();
        const std::vector<libMesh::Point>& xyz = fe->get_xyz();
        const unsigned int n_qp = (unsigned int)fe->get_qpoints().size();

        // Initialize some variables
        Real rho = 0.0;
        Real drho = 0.0;
        Real Vc = 0.0;
        Real dVc = 0.0;
        const MAST::FieldFunction<Real>* Vc_f;

        // Create a unity field function that evaluate to 1 over the entire domain, used for 3D elements
        MAST::Parameter one("one", 1.0);
        MAST::ConstantFieldFunction one_f("one", one);

        // Create a field function to calculate the area of a 1D rectangular element (old MAST)
        class Area1DFieldFunction: public MAST::FieldFunction<Real> {
        public:
            Area1DFieldFunction(const MAST::FieldFunction<Real>& hy, const MAST::FieldFunction<Real>& hz):
            MAST::FieldFunction<Real>("A"), _hy(hy), _hz(hz) {}
            virtual void operator()(const libMesh::Point& p, const Real t, Real& v) const {
                Real hy;
                Real hz;
                _hy(p, t, hy);
                _hz(p, t, hz);
                v = hy * hz;
            }
            virtual void derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t, Real& v) const {
                Real hy, hz, dhy, dhz;
                if (&f == &_hy) {
                    _hz(p, t, hz);
                    v = hz;
                }
                else if (&f == &_hz) {
                    _hy(p, t, hy);
                    v = hy;
                }
                else {
                    _hy(p, t, hy);
                    _hy.derivative(f, p, t, dhy);
                    _hz(p, t, hz);
                    _hz.derivative(f, p, t, dhz);
                    v = dhy * hz + hy * dhz;
                }
            }
        protected:
            const MAST::FieldFunction<Real>& _hy;
            const MAST::FieldFunction<Real>& _hz;
        };
        Area1DFieldFunction A1d_f(section_property.get<MAST::FieldFunction<Real>>("hy"),
                                  section_property.get<MAST::FieldFunction<Real>>("hz"));
        // TODO: Once the 1D cross-section library is implemented in this branch, the class and field function above will no longer be needed.
        //   this is only used for "old" MAST where only rectangular cross sections were supported.

        // To get the total volume for integration, we need to multiply by area for 1D elements and by thickness for 2D elements
        if (section_property.dim() == 1) {
            // Get the area for 1D elements since volume = JxW * A in this case
            // Vc_f = &(section_property.get<MAST::FieldFunction<Real>>("A"));
            Vc_f = &A1d_f;
        }
        else if (section_property.dim() == 2) {
            // Get the thickness for 2D elements since volume = JxW * h in this case
            Vc_f = &(section_property.get<MAST::FieldFunction<Real>>("h"));
        }
        else if (section_property.dim() == 3) {
            // Get unity (1) for 3D element since volume = JxW in this case
            Vc_f = &one_f;
        }

        Real time = 0.0;
        // Iterate over the integration points to perform the numerical integration
        for (unsigned int ind_qp=0; ind_qp<n_qp; ind_qp++) {
            // Get the current time of the system
            time = _system->system().time;
            // Get the density at this integration point
            rho_f(xyz[ind_qp], time, rho);
            rho_f.derivative(f, xyz[ind_qp], time, drho); 
            // Get the extra volume component at this integration point
            (*Vc_f)(xyz[ind_qp], time, Vc);
            Vc_f->derivative(f, xyz[ind_qp], time, dVc);
            Real dxyz = 0.0;  // FIXME: If the parameter is nodal coordinates (shape sensitivity) then this may be nonzero
            _dfirst_moment_of_inertia_yz_dp += (drho * xyz[ind_qp](0) * Vc + rho * dxyz * Vc + rho * xyz[ind_qp](0) * dVc) * JxW[ind_qp];
            _dfirst_moment_of_inertia_xz_dp += (drho * xyz[ind_qp](1) * Vc + rho * dxyz * Vc + rho * xyz[ind_qp](1) * dVc) * JxW[ind_qp];
            _dfirst_moment_of_inertia_xy_dp += (drho * xyz[ind_qp](2) * Vc + rho * dxyz * Vc + rho * xyz[ind_qp](2) * dVc) * JxW[ind_qp];
            _dmass_dp += (drho * Vc + rho * dVc) * JxW[ind_qp];
        }
    }
}


void
MAST::CenterOfMassOutput::
evaluate_topology_sensitivity(const MAST::FunctionBase &f) {
    
    // the primal data should have been calculated
    libmesh_assert(_physics_elem);
    libmesh_assert(f.is_topology_parameter());

    libmesh_error_msg("Must be implemented in derived class; " << __PRETTY_FUNCTION__ << " in " << __FILE__ << " at line number " << __LINE__);
    
    // if (this->if_evaluate_for_element(_physics_elem->elem())) {
        
    //     std::pair<const MAST::FieldFunction<RealVectorX>*, unsigned int>
    //     val = this->get_elem_boundary_velocity_data();
        
    //     if (val.first) {

    //         MAST::StructuralElementBase& e =
    //         dynamic_cast<MAST::StructuralElementBase&>(*_physics_elem);
            
    //         RealVectorX
    //         vec   = RealVectorX::Zero(e.sol().size());
            
    //         RealMatrixX
    //         dummy = RealMatrixX::Zero(vec.size(), vec.size());
            
    //         e.volume_external_residual_boundary_velocity(f,
    //                                                      val.second,
    //                                                      *val.first,
    //                                                      _discipline->volume_loads(),
    //                                                      false,
    //                                                      vec,
    //                                                      dummy);
            
    //         // compute the contribution of this element to center_of_mass
    //         _dcenter_of_mass_dp -= vec.dot(e.sol());
    //     }
    // }
}



void
MAST::CenterOfMassOutput::
evaluate_topology_sensitivity(const MAST::FunctionBase &f,
                              const MAST::FieldFunction<RealVectorX> &vel) {
    
    // the primal data should have been calculated
    libmesh_assert(_physics_elem);
    libmesh_assert(f.is_topology_parameter());

    libmesh_error_msg("Must be implemented in derived class; " << __PRETTY_FUNCTION__ << " in " << __FILE__ << " at line number " << __LINE__);

    // const MAST::LevelSetIntersectedElem
    // &elem = dynamic_cast<const MAST::LevelSetIntersectedElem&>(_physics_elem->elem());
    
    // // sensitivity only exists at the boundary. So, we proceed with calculation
    // // only if this element has an intersection in the interior, or with a side.
    
    // if (this->if_evaluate_for_element(elem)    &&
    //     elem.if_elem_has_level_set_boundary()       &&
    //     elem.if_subelem_has_side_on_level_set_boundary()) {
        
    //     MAST::StructuralElementBase& e =
    //     dynamic_cast<MAST::StructuralElementBase&>(*_physics_elem);
        
    //     RealVectorX
    //     vec   = RealVectorX::Zero(e.sol().size());
        
    //     RealMatrixX
    //     dummy = RealMatrixX::Zero(vec.size(), vec.size());
        
    //     e.volume_external_residual_boundary_velocity(f,
    //                                                  elem.get_subelem_side_on_level_set_boundary(),
    //                                                  vel,
    //                                                  _discipline->volume_loads(),
    //                                                  false,
    //                                                  vec,
    //                                                  dummy);
        
    //     // compute the contribution of this element to center_of_mass
    //     _dcenter_of_mass_dp -= vec.dot(e.sol());
    // }
}



RealVectorX
MAST::CenterOfMassOutput::output_vector_total() {
    
    RealVectorX val = _center_of_mass;
    Real first_moment_of_inertia_yz = _first_moment_of_inertia_yz;
    Real first_moment_of_inertia_xz = _first_moment_of_inertia_xz;
    Real first_moment_of_inertia_xy = _first_moment_of_inertia_xy;
    Real mass = _mass;
    
    if (!_skip_comm_sum) {
        _system->system().comm().sum(first_moment_of_inertia_yz);
        _system->system().comm().sum(first_moment_of_inertia_xz);
        _system->system().comm().sum(first_moment_of_inertia_xy);
        _system->system().comm().sum(mass);
    }

    val.coeffRef(0,0) = first_moment_of_inertia_yz / mass;
    val.coeffRef(1,0) = first_moment_of_inertia_xz / mass;
    val.coeffRef(2,0) = first_moment_of_inertia_xy / mass;
    
    return val;
}



RealVectorX
MAST::CenterOfMassOutput::output_vector_sensitivity_total(const MAST::FunctionBase& p) {

    Real first_moment_of_inertia_yz = _first_moment_of_inertia_yz;
    Real first_moment_of_inertia_xz = _first_moment_of_inertia_xz;
    Real first_moment_of_inertia_xy = _first_moment_of_inertia_xy;
    Real mass = _mass;

    RealVectorX val = _dcenter_of_mass_dp;
    Real dfirst_moment_of_inertia_yz_dp = _dfirst_moment_of_inertia_yz_dp;
    Real dfirst_moment_of_inertia_xz_dp = _dfirst_moment_of_inertia_xz_dp;
    Real dfirst_moment_of_inertia_xy_dp = _dfirst_moment_of_inertia_xy_dp;
    Real dmass_dp = _dmass_dp;
    
    if (!_skip_comm_sum) {
        _system->system().comm().sum(first_moment_of_inertia_yz);
        _system->system().comm().sum(first_moment_of_inertia_xz);
        _system->system().comm().sum(first_moment_of_inertia_xy);
        _system->system().comm().sum(mass);

        _system->system().comm().sum(dfirst_moment_of_inertia_yz_dp);
        _system->system().comm().sum(dfirst_moment_of_inertia_xz_dp);
        _system->system().comm().sum(dfirst_moment_of_inertia_xy_dp);
        _system->system().comm().sum(dmass_dp);
    }

    val.coeffRef(0,0) = (dfirst_moment_of_inertia_yz_dp - first_moment_of_inertia_yz * dmass_dp / mass) / mass;
    val.coeffRef(1,0) = (dfirst_moment_of_inertia_xz_dp - first_moment_of_inertia_xz * dmass_dp / mass) / mass;
    val.coeffRef(2,0) = (dfirst_moment_of_inertia_xy_dp - first_moment_of_inertia_xy * dmass_dp / mass) / mass;

    return val;
}



void
MAST::CenterOfMassOutput::output_derivative_for_elem(RealVectorX& dq_dX) {
    
    // make sure that this has not been initialized ana calculated for all elems
    libmesh_assert(_physics_elem);

    libmesh_error_msg("Must be implemented in derived class; " << __PRETTY_FUNCTION__ << " in " << __FILE__ << " at line number " << __LINE__);

    if (this->if_evaluate_for_element(_physics_elem->elem())) {
        dq_dX.setZero();
        // The center of mass is a function of displacements so I need to figure out the correct equation for this
    }
}


void
MAST::CenterOfMassOutput::
set_elem_data(unsigned int dim,
              const libMesh::Elem& ref_elem,
              MAST::GeomElem& elem) const {
    
    libmesh_assert(!_physics_elem);
}


void
MAST::CenterOfMassOutput::init(const MAST::GeomElem& elem) {
    
    libmesh_assert(!_physics_elem);
    libmesh_assert(_assembly);
    libmesh_assert(_system);
    const MAST::ElementPropertyCardBase& p =
    dynamic_cast<const MAST::ElementPropertyCardBase&>(_discipline->get_property_card(elem));
    
    _physics_elem =
    MAST::build_structural_element(*_system, elem, p).release(); // TODO: Should this be a structural element?
}


