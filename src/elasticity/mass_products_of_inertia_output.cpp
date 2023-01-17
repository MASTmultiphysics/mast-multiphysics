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
#include "elasticity/mass_products_of_inertia_output.h"
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


MAST::MassProductsOfInertiaOutput::MassProductsOfInertiaOutput():
MAST::OutputAssemblyElemOperations() {
    _mass_products_of_inertia = RealVectorX::Zero(3,1);
    _dmass_products_of_inertia_dp = RealVectorX::Zero(3,1);
    _product_of_inertia_xy = 0.0;
    _product_of_inertia_yz = 0.0;
    _product_of_inertia_zx = 0.0;
    _dproduct_of_inertia_xy_dp = 0.0;
    _dproduct_of_inertia_yz_dp = 0.0;
    _dproduct_of_inertia_zx_dp = 0.0;
}


MAST::MassProductsOfInertiaOutput::~MassProductsOfInertiaOutput() {

}


void
MAST::MassProductsOfInertiaOutput::zero_for_analysis() {
    _mass_products_of_inertia.setZero();
    _dmass_products_of_inertia_dp.setZero();
    _product_of_inertia_xy = 0.0;
    _product_of_inertia_yz = 0.0;
    _product_of_inertia_zx = 0.0;
    _dproduct_of_inertia_xy_dp = 0.0;
    _dproduct_of_inertia_yz_dp = 0.0;
    _dproduct_of_inertia_zx_dp = 0.0;
}


void
MAST::MassProductsOfInertiaOutput::zero_for_sensitivity() {
    _mass_products_of_inertia.setZero();
    _dmass_products_of_inertia_dp.setZero();
    _product_of_inertia_xy = 0.0;
    _product_of_inertia_yz = 0.0;
    _product_of_inertia_zx = 0.0;
    _dproduct_of_inertia_xy_dp = 0.0;
    _dproduct_of_inertia_yz_dp = 0.0;
    _dproduct_of_inertia_zx_dp = 0.0;
}


void
MAST::MassProductsOfInertiaOutput::evaluate() {
    
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
        // TODO: Once the 1D cross-section library is implemented in this branch, the class and field function above will no longer be needed.
        //   this is only used for "old" MAST where only rectangular cross sections were supported.

        // To get the total volume for integration, we need to multiply by area for 1D elements and by thickness for 2D elements
        std::unique_ptr<Area1DFieldFunction> A1d_f;
        const MAST::FieldFunction<Real>* Vc_f;
        if (section_property.dim() == 1) {
            A1d_f = std::unique_ptr<Area1DFieldFunction>(new Area1DFieldFunction(section_property.get<MAST::FieldFunction<Real>>("hy"),
                                                                                 section_property.get<MAST::FieldFunction<Real>>("hz")));
            // Get the area for 1D elements since volume = JxW * A in this case
            // Vc_f = &(section_property.get<MAST::FieldFunction<Real>>("A"));
            Vc_f = A1d_f.get();
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
        Real x, y, z;
        // Iterate over the integration points to perform the numerical integration
        for (unsigned int ind_qp=0; ind_qp<n_qp; ind_qp++) {
            // Get the current time of the system
            time = _system->system().time;

            // Get the density at this integration point
            rho_f(xyz[ind_qp], time, rho);
            // Get the extra volume component at this integration point
            (*Vc_f)(xyz[ind_qp], time, Vc);

            mass_elem = rho * Vc * JxW[ind_qp]; // A value used in multiple calculations below
            x = xyz[ind_qp](0);
            y = xyz[ind_qp](1);
            z = xyz[ind_qp](2);

            _mass_products_of_inertia.coeffRef(0,0) += mass_elem * (x * y);
            _mass_products_of_inertia.coeffRef(1,0) += mass_elem * (y * z);
            _mass_products_of_inertia.coeffRef(2,0) += mass_elem * (x * z);

            _product_of_inertia_xy += mass_elem * (x * y);
            _product_of_inertia_yz += mass_elem * (y * z);
            _product_of_inertia_zx += mass_elem * (x * z);
        }
    }
    // FIXME: Product moments of inertia will change with displacements. Need to account for that in these calculations.
}


void
MAST::MassProductsOfInertiaOutput::evaluate_sensitivity(const MAST::FunctionBase &f) {
    
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
        // TODO: Once the 1D cross-section library is implemented in this branch, the class and field function above will no longer be needed.
        //   this is only used for "old" MAST where only rectangular cross sections were supported.

        // To get the total volume for integration, we need to multiply by area for 1D elements and by thickness for 2D elements
        std::unique_ptr<Area1DFieldFunction> A1d_f;
        const MAST::FieldFunction<Real>* Vc_f;
        if (section_property.dim() == 1) {
            A1d_f = std::unique_ptr<Area1DFieldFunction>(new Area1DFieldFunction(section_property.get<MAST::FieldFunction<Real>>("hy"),
                                                                                 section_property.get<MAST::FieldFunction<Real>>("hz")));
            // Get the area for 1D elements since volume = JxW * A in this case
            // Vc_f = &(section_property.get<MAST::FieldFunction<Real>>("A"));
            Vc_f = A1d_f.get();
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
        Real x, y, z, dx, dy, dz;
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
            x = xyz[ind_qp](0);
            y = xyz[ind_qp](1);
            z = xyz[ind_qp](2);
            dx = 0;  // FIXME: If the parameter is nodal coordinates (shape sensitivity) then this may be nonzero
            dy = 0;  // FIXME: If the parameter is nodal coordinates (shape sensitivity) then this may be nonzero
            dz = 0;  // FIXME: If the parameter is nodal coordinates (shape sensitivity) then this may be nonzero

            _dmass_products_of_inertia_dp.coeffRef(0,0) += (dx*y*rho*Vc + x*dy*rho*Vc + x*y*drho*Vc + x*y*rho*dVc) * JxW[ind_qp];
            _dmass_products_of_inertia_dp.coeffRef(1,0) += (dy*z*rho*Vc + y*dz*rho*Vc + y*z*drho*Vc + y*z*rho*dVc) * JxW[ind_qp];
            _dmass_products_of_inertia_dp.coeffRef(2,0) += (dz*x*rho*Vc + z*dx*rho*Vc + z*x*drho*Vc + z*x*rho*dVc) * JxW[ind_qp];

            _dproduct_of_inertia_xy_dp += (dx*y*rho*Vc + x*dy*rho*Vc + x*y*drho*Vc + x*y*rho*dVc) * JxW[ind_qp];
            _dproduct_of_inertia_yz_dp += (dy*z*rho*Vc + y*dz*rho*Vc + y*z*drho*Vc + y*z*rho*dVc) * JxW[ind_qp];
            _dproduct_of_inertia_zx_dp += (dz*x*rho*Vc + z*dx*rho*Vc + z*x*drho*Vc + z*x*rho*dVc) * JxW[ind_qp];
        }
    }
}


void
MAST::MassProductsOfInertiaOutput::
evaluate_topology_sensitivity(const MAST::FunctionBase &f) {
    
    // the primal data should have been calculated
    libmesh_assert(_physics_elem);
    libmesh_assert(f.is_topology_parameter());

    libmesh_error_msg("Must be implemented in derived class; " << __PRETTY_FUNCTION__ << " in " << __FILE__ << " at line number " << __LINE__);
}



void
MAST::MassProductsOfInertiaOutput::
evaluate_topology_sensitivity(const MAST::FunctionBase &f,
                              const MAST::FieldFunction<RealVectorX> &vel) {
    
    // the primal data should have been calculated
    libmesh_assert(_physics_elem);
    libmesh_assert(f.is_topology_parameter());

    libmesh_error_msg("Must be implemented in derived class; " << __PRETTY_FUNCTION__ << " in " << __FILE__ << " at line number " << __LINE__);
}



RealVectorX
MAST::MassProductsOfInertiaOutput::output_vector_total() {
    
    RealVectorX val = _mass_products_of_inertia;
    Real product_of_inertia_xy = _product_of_inertia_xy;
    Real product_of_inertia_yz = _product_of_inertia_yz;
    Real product_of_inertia_zx = _product_of_inertia_zx;
    
    if (!_skip_comm_sum) {
        _system->system().comm().sum(product_of_inertia_xy);
        _system->system().comm().sum(product_of_inertia_yz);
        _system->system().comm().sum(product_of_inertia_zx);
    }

    val.coeffRef(0,0) = product_of_inertia_xy;
    val.coeffRef(1,0) = product_of_inertia_yz;
    val.coeffRef(2,0) = product_of_inertia_zx;
    
    return val;
}



RealVectorX
MAST::MassProductsOfInertiaOutput::output_vector_sensitivity_total(const MAST::FunctionBase& p) {

    Real product_of_inertia_xy = _product_of_inertia_xy;
    Real product_of_inertia_yz = _product_of_inertia_yz;
    Real product_of_inertia_zx = _product_of_inertia_zx;

    RealVectorX val = _dmass_products_of_inertia_dp;
    Real dproduct_of_inertia_xy_dp = _dproduct_of_inertia_xy_dp;
    Real dproduct_of_inertia_yz_dp = _dproduct_of_inertia_yz_dp;
    Real dproduct_of_inertia_zx_dp = _dproduct_of_inertia_zx_dp;
    
    if (!_skip_comm_sum) {
        _system->system().comm().sum(product_of_inertia_xy);
        _system->system().comm().sum(product_of_inertia_yz);
        _system->system().comm().sum(product_of_inertia_zx);

        _system->system().comm().sum(dproduct_of_inertia_xy_dp);
        _system->system().comm().sum(dproduct_of_inertia_yz_dp);
        _system->system().comm().sum(dproduct_of_inertia_zx_dp);
    }

    val.coeffRef(0,0) = dproduct_of_inertia_xy_dp;
    val.coeffRef(1,0) = dproduct_of_inertia_yz_dp;
    val.coeffRef(2,0) = dproduct_of_inertia_zx_dp;

    return val;
}



void
MAST::MassProductsOfInertiaOutput::output_derivative_for_elem(RealVectorX& dq_dX) {
    
    // make sure that this has not been initialized ana calculated for all elems
    libmesh_assert(_physics_elem);

    libmesh_error_msg("Must be implemented in derived class; " << __PRETTY_FUNCTION__ << " in " << __FILE__ << " at line number " << __LINE__);

    if (this->if_evaluate_for_element(_physics_elem->elem())) {
        dq_dX.setZero();
        // The center of mass is a function of displacements so I need to figure out the correct equation for this
    }
}


void
MAST::MassProductsOfInertiaOutput::
set_elem_data(unsigned int dim,
              const libMesh::Elem& ref_elem,
              MAST::GeomElem& elem) const {
    
    libmesh_assert(!_physics_elem);
}


void
MAST::MassProductsOfInertiaOutput::init(const MAST::GeomElem& elem) {
    
    libmesh_assert(!_physics_elem);
    libmesh_assert(_assembly);
    libmesh_assert(_system);
    const MAST::ElementPropertyCardBase& p =
    dynamic_cast<const MAST::ElementPropertyCardBase&>(_discipline->get_property_card(elem));
    
    _physics_elem =
    MAST::build_structural_element(*_system, elem, p).release(); // TODO: Should this be a structural element?
}


