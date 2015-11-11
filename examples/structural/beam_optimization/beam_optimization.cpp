/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2015  Manav Bhatia
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

// C++ includes
#include <iostream>

// MAST includes
#include "examples/structural/beam_optimization/beam_optimization.h"
#include "driver/driver_base.h"

// libMesh includes
#include "libmesh/parallel_mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"



MAST::MultilinearInterpolation::MultilinearInterpolation
(const std::string& nm,
 std::map<Real, MAST::FieldFunction<Real>*>& values):
MAST::FieldFunction<Real>(nm),
_values(values) {
    
    // make sure that the size of the provided values is finite
    libmesh_assert(values.size() > 0);
    
    std::map<Real, MAST::FieldFunction<Real>*>::iterator
    it = values.begin(), end = values.end();
    
    // tell the function that it is dependent on the provided functions
    for ( ; it != end; it++)
        _functions.insert(it->second->master());
}

MAST::MultilinearInterpolation::MultilinearInterpolation
(const MAST::MultilinearInterpolation& o):
MAST::FieldFunction<Real>(o),
_values(o._values) {
    std::map<Real, MAST::FieldFunction<Real>*>::iterator
    it = _values.begin(), end = _values.end();
    
    // tell the function that it is dependent on the provided functions
    for ( ; it != end; it++)
        _functions.insert(it->second->master());
}


std::auto_ptr<MAST::FieldFunction<Real> >
MAST::MultilinearInterpolation::clone() const {
    
    return std::auto_ptr<MAST::FieldFunction<Real> >
    (new MAST::MultilinearInterpolation(*this));
}

MAST::MultilinearInterpolation::~MultilinearInterpolation() {
    
}


void
MAST::MultilinearInterpolation::operator() (const libMesh::Point& p,
                                            Real t,
                                            Real& v) const {
    
    //
    // the following is used for calculation of the return value
    //   f(x) is defined for x for each x0 < x < x1
    //   if   x <= x0,      f(x) = f(x0)
    //   if   x0 < x < x1,  f(x) is interpolated
    //   if   x >= x1,      f(x) = f(x1)
    //
    
    std::map<Real, MAST::FieldFunction<Real>*>::const_iterator
    it1, it2;
    std::map<Real, MAST::FieldFunction<Real>*>::const_reverse_iterator
    rit = _values.rbegin();
    it1  = _values.begin();
    
    // check the lower bound
    if (p(0) <=  it1->first) {
        (*it1->second)(p, t, v);
    }
    // check the upper bound
    else if (p(0) >=  rit->first) {
        (*rit->second)(p, t, v);
    }
    else {
        // if it gets here, the ordinate is in between the provided range
        it2 = _values.lower_bound(p(0));
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
              (p(0) - it1->first)/(it2->first - it1->first) *
              (f1-f0));
    }
}

void
MAST::MultilinearInterpolation::derivative(const MAST::DerivativeType d,
                                           const MAST::FunctionBase& f,
                                           const libMesh::Point& p,
                                           Real t,
                                           Real& v) const {
    
    //
    // the following is used for calculation of the return value
    //   f(x) is defined for x for each x0 < x < x1
    //   if   x <= x0,      f(x) = f(x0)
    //   if   x0 < x < x1,  f(x) is interpolated
    //   if   x >= x1,      f(x) = f(x1)
    //
    
    std::map<Real, MAST::FieldFunction<Real>*>::const_iterator
    it1, it2;
    std::map<Real, MAST::FieldFunction<Real>*>::const_reverse_iterator
    rit = _values.rbegin();
    it1  = _values.begin();
    
    // check the lower bound
    if (p(0) <=  it1->first) {
        (*it1->second)(p, t, v);
    }
    // check the upper bound
    else if (p(0) >=  rit->first) {
        (*rit->second)(p, t, v);
    }
    else {
        // if it gets here, the ordinate is in between the provided range
        it2 = _values.lower_bound(p(0));
        // this cannot be the first element of the map
        libmesh_assert(it2 != _values.begin());
        // it2 provides the upper bound. The lower bound is provided by the
        // preceding iterator
        it1 = it2--;
        Real f0 = 0., f1 = 0.;
        it1->second->derivative(d, f, p, t, f0);
        it2->second->derivative(d, f, p, t, f1);
        // now interpolate
        v =  (f0 +
              (p(0) - it1->first)/(it2->first - it1->first) *
              (f1-f0));
    }
}




MAST::BeamOffset::BeamOffset(const std::string& nm,
                             MAST::FieldFunction<Real> *thickness):
MAST::FieldFunction<Real>(nm),
_dim(thickness) {
    
    _functions.insert(thickness->master());
}

MAST::BeamOffset::BeamOffset(const MAST::BeamOffset& o):
MAST::FieldFunction<Real>(o),
_dim(o._dim->clone().release()) {
    
    _functions.insert(_dim->master());
}


std::auto_ptr<MAST::FieldFunction<Real> >
MAST::BeamOffset::clone() const {
    
    return std::auto_ptr<MAST::FieldFunction<Real> >
    (new MAST::BeamOffset(*this));
}


MAST::BeamOffset::~BeamOffset() {
    delete _dim;
}


void
MAST::BeamOffset::operator() (const libMesh::Point& p,
                              Real t,
                              Real& v) const {
    
    (*_dim)(p, t, v);
    v *= 0.5;
}


void
MAST::BeamOffset::derivative(const MAST::DerivativeType d,
                             const MAST::FunctionBase& f,
                             const libMesh::Point& p,
                             Real t,
                             Real& v) const {
    _dim->derivative(d, f, p, t, v);
    v *= 0.5;
}




MAST::Weight::Weight(MAST::PhysicsDisciplineBase& discipline):
MAST::FieldFunction<Real>("Weight"),
_discipline(discipline)
{ }



MAST::Weight::Weight(const MAST::Weight& w):
MAST::FieldFunction<Real>(w),
_discipline(w._discipline)
{ }



std::auto_ptr<MAST::FieldFunction<Real> >
MAST::Weight::clone() const {
    
    return std::auto_ptr<MAST::FieldFunction<Real> >
    (new MAST::Weight(*this));
}



MAST::Weight::~Weight() { }



void
MAST::Weight::operator() (const libMesh::Point& p,
                          Real t,
                          Real& v) const {
    
    // get a reference to the mesh
    const libMesh::MeshBase&
    mesh    = _discipline.get_equation_systems().get_mesh();
    libMesh::MeshBase::const_element_iterator
    eit     = mesh.active_local_elements_begin(),
    eend    = mesh.active_local_elements_end();
    
    Real h, rho, x0, x1, dx;
    v = 0.;
    
    libMesh::Point elem_p;
    const unsigned int n_sec = 3; // number of quadrature divs
    
    for ( ; eit != eend; eit++ ) {
        
        // get a pointer to the element and then as the discipline
        // for the element property
        const libMesh::Elem* e = *eit;
        
        const MAST::ElementPropertyCardBase& prop =
        _discipline.get_property_card(*e);
        
        // assuming that the element is one-dimensional, we need
        // its section area value
        
        // before that, convert the property to a 1D section property
        // card
        const MAST::Solid1DSectionElementPropertyCard& prop1d =
        dynamic_cast<const MAST::Solid1DSectionElementPropertyCard&>(prop);
        
        // get a reference to the section area
        const MAST::FieldFunction<Real>& area = prop1d.A();
        
        // get a reference to the density variable
        const MAST::MaterialPropertyCardBase& mat = prop.get_material();
        const MAST::FieldFunction<Real> &rhof =
        mat.get<MAST::FieldFunction<Real> >("rho");
        
        
        // for each element iterate over the length and calculate the
        // weight from the section area and section density
        // use three point trapezoidal rule to calculate the integral
        x0 = e->point(0)(0);
        x1 = e->point(1)(0);
        dx = (x1-x0)/n_sec;
        for (unsigned int i=0; i<n_sec; i++) {
            elem_p(0) = x0 + dx*(i+0.5);
            area(elem_p, 0., h);
            rhof(elem_p, 0., rho);
            v += h * rho * dx;
        }
        
    }
}





void
MAST::Weight::derivative(const MAST::DerivativeType d,
                         const MAST::FunctionBase& f,
                         const libMesh::Point& p,
                         Real t,
                         Real& v) const {
    
    // get a reference to the mesh
    const libMesh::MeshBase&
    mesh    = _discipline.get_equation_systems().get_mesh();
    libMesh::MeshBase::const_element_iterator
    eit     = mesh.active_local_elements_begin(),
    eend    = mesh.active_local_elements_end();
    
    Real h, rho, dh, drho, x0, x1, dx;
    v = 0.;
    
    libMesh::Point elem_p;
    const unsigned int n_sec = 3; // number of quadrature divs
    
    for ( ; eit != eend; eit++ ) {
        
        // get a pointer to the element and then as the discipline
        // for the element property
        const libMesh::Elem* e = *eit;
        
        const MAST::ElementPropertyCardBase& prop =
        _discipline.get_property_card(*e);
        
        // assuming that the element is one-dimensional, we need
        // its section area value
        
        // before that, convert the property to a 1D section property
        // card
        const MAST::Solid1DSectionElementPropertyCard& prop1d =
        dynamic_cast<const MAST::Solid1DSectionElementPropertyCard&>(prop);
        
        // get a reference to the section area
        const MAST::FieldFunction<Real>&
        area = prop1d.A();
        
        // get a reference to the density variable
        const MAST::MaterialPropertyCardBase& mat =
        prop.get_material();
        const MAST::FieldFunction<Real> &rhof =
        mat.get<MAST::FieldFunction<Real> >("rho");
        
        
        // for each element iterate over the length and calculate the
        // weight from the section area and section density
        // use three point trapezoidal rule to calculate the integral
        x0 = e->point(0)(0);
        x1 = e->point(1)(0);
        dx = (x1-x0)/n_sec;
        for (unsigned int i=0; i<n_sec; i++) {
            elem_p(0) = x0 + dx*(i+0.5);
            area(elem_p, 0., h);
            area.derivative(d, f, elem_p, 0., dh);
            rhof(elem_p, 0., rho);
            rhof.derivative(d, f, elem_p, 0., drho);
            v += (dh * rho + h * drho) * dx;
        }
        
    }
}




MAST::BeamBendingSizingOptimization::
BeamBendingSizingOptimization(libMesh::LibMeshInit& init,
                              std::ostream& output):
MAST::FunctionEvaluation(output),
_libmesh_init(init),
_n_elems(0),
_n_stations(0),
_disp_0(0.),
_vf_0(0.),
_weight(NULL),
_mesh(NULL),
_eq_systems(NULL),
_static_system(NULL),
_structural_discipline(NULL),
_structural_sys(NULL),
_dirichlet(NULL),
_flux_load(NULL) {
    
    // call the initialization routine to setup the data structures
    _init();
}




MAST::BeamBendingSizingOptimization::~BeamBendingSizingOptimization() {
    
    delete   _weight;
    delete   _eq_systems;
    delete   _mesh;
    
    
    // delete the h_y station values
    std::map<Real, MAST::FieldFunction<Real>*>::iterator
    it  = _h_y_station_vals.begin(),
    end = _h_y_station_vals.end();
    for (; it != end; it++)
        delete it->second;
    
    for (unsigned int i=0; i<_n_vars; i++) {
        delete _disp_function_sens[i];
    }
}



void
MAST::BeamBendingSizingOptimization::init_dvar(std::vector<Real>& x,
                                               std::vector<Real>& xmin,
                                               std::vector<Real>& xmax) {
    // one DV for each element
    x       = _dv_init;
    xmin    = _dv_low;
    xmax.resize(_n_vars);
    std::fill(xmax.begin(), xmax.end(), 1.);
}



void
MAST::BeamBendingSizingOptimization::evaluate(const std::vector<Real>& dvars,
                                              Real& obj,
                                              bool eval_obj_grad,
                                              std::vector<Real>& obj_grad,
                                              std::vector<Real>& fvals,
                                              std::vector<bool>& eval_grads,
                                              std::vector<Real>& grads) {
    
    
    libmesh_assert_equal_to(dvars.size(), _n_vars);
    
    // set the parameter values equal to the DV value
    for (unsigned int i=0; i<_n_vars; i++)
        *_parameters[i] = dvars[i]*_dv_scaling[i];
    
    // DO NOT zero out the gradient vector, since GCMMA needs it for the
    // subproblem solution
    //std::fill(obj_grad.begin(), obj_grad.end(), 0.);
    //std::fill(grads.begin(), grads.end(), 0.);
    
    libMesh::Point pt; // dummy point object
    
    libMesh::out << "New Eval" << std::endl;
    
    // the optimization problem is defined as
    // min weight, subject to constraints on displacement and stresses
    Real wt = 0., vf = 0., disp = 0.;
    
    
    // calculate weight
    (*_weight)(pt, 0., wt);
    
    
    // first zero the solution and init the Euler variables to undisturbed solution
    _static_system->solution->zero();
    
    // now solve for this load step
    _static_system->solve();
    {
        std::set<std::string> names;
        names.insert("StaticStructuralSystem");
        libMesh::ExodusII_IO(*_mesh).write_equation_systems("str.exo",
                                                            *_eq_systems,
                                                            &names);
    }
    
    
    // now get the displacement constraint
    pt(0) = 3.;
    DenseRealVector disp_vec;
    //(*_disp_function)(pt, 0., disp_vec);
    // reference displacement value
    // w < w0 => w/w0 < 1. => w/w0 - 1. < 0
    //disp = disp_vec(0);
    
    std::vector<Real> grad_vals;
    
    // set the function and objective values
    // flutter objective
    // vf > v0 => vf/v0 > 1 => 1-vf/v0 < 0
    //
    obj = wt;
    fvals[0] = disp/_disp_0-1.;
    fvals[1] =     1.-vf/_vf_0;
    Real w_sens = 0.;
    
    // set gradient of weight
    for (unsigned int i=0; i<_n_vars; i++) {
        _weight->derivative(MAST::PARTIAL_DERIVATIVE,
                            *_parameter_functions[i],
                            pt,
                            0.,
                            w_sens);
        obj_grad[i] = w_sens*_dv_scaling[i];
    }
    
    
    // evaluate sensitivity if any of the function sensitivity is required
    bool if_sens = (false || eval_obj_grad);
    for (unsigned int i=0; i<eval_grads.size(); i++)
        if_sens = (if_sens || eval_grads[i]);
    if_sens = false;
    if (if_sens) {
        // grad_k = dfi/dxj  ,  where k = j*NFunc + i
        
        // each design variable requires its own set of iterations to converge
        // the sensitivity
        // the design variables are iterated upon in the reverse order since
        // we are doing this one DV at a time, and if libMesh::ParameterVector
        // has size one, libMesh will put the design sensitivity vector in
        // System::get_sensitivity_solution(0).
        // add the sensitivity solution for param num 0
        _static_system->add_sensitivity_solution(0);
        
        // get the displacement gradient
        pt(0) = 3.0;
        for (unsigned int j=0; j<_n_vars; j++) {
            disp_vec.zero();
            //(*_disp_function_sens[j])(pt, 0., disp_vec);
            //grads[j*_n_ineq] = disp_vec(0)/_disp_0*_dv_scaling[j];
        }
    }
    
    // write the evaluation output
    this->output(0, dvars, obj, fvals, false);
}



void
MAST::BeamBendingSizingOptimization::_init() {
    
    // use input.in as the input file
    GetPot infile("input.in");
    
    _mesh                  = new libMesh::ParallelMesh(_libmesh_init.comm());
    
    // design data
    _n_elems               = infile("n_elems", 50);
    _n_stations            = infile("n_stations", 2);
    _n_vars                = _n_stations; // for thickness variable
    
    _n_eq                  = 0;
    _n_ineq                = (1 +          // +1 for the displacement
                              _n_elems);   // +1 for element stress
    _max_iters             = 10000;
    
    // initialize the dv vector data
    Real
    th_l                   = infile("thickness_lower", 0.001),
    th_u                   = infile("thickness_upper", 0.2),
    th                     = infile("thickness", 0.01),
    
    // displacement constraint
    _disp_0                = infile( "displacement_0", 0.);
    
    _dv_init.resize    (_n_vars);
    _dv_scaling.resize (_n_vars);
    _dv_low.resize     (_n_vars);
    
    // design variables for the thickness values
    for (unsigned int i=0; i<_n_vars; i++) {
        
        _dv_init[i]    =  infile("dv_init", th/th_u, i);
        _dv_low[i]     = th_l/th_u;
        _dv_scaling[i] =      th_u;
    }
    
    
    // now initialize the mesh
    libMesh::MeshTools::Generation::build_line(*_mesh, _n_elems);
    
    // Print information about the mesh to the screen.
    _mesh->print_info();
    
    // Create an equation systems object.
    _eq_systems            = new libMesh::EquationSystems(*_mesh);
    
    // Declare the system
    _static_system         =
    &(_eq_systems->add_system<libMesh::NonlinearImplicitSystem>
      ("StaticStructuralSystem"));
    
    libMesh::FEType fetype (libMesh::FIRST, libMesh::LAGRANGE);
    
    _structural_discipline.reset(new MAST::StructuralDiscipline(*_eq_systems));
    _structural_sys.reset(new MAST::StructuralSystemInitialization(*_static_system,
                                                                   _static_system->name(),
                                                                   fetype));
    
    
    // set the pressure loads
    _pressure.reset(new MAST::Parameter("p", infile("pressure", 2.)));
    _pressure_fn.reset(new MAST::ConstantFieldFunction("pressure", *_pressure));
    _flux_load.reset(new MAST::BoundaryConditionBase(MAST::SURFACE_PRESSURE));
    _flux_load->add(*_pressure_fn);
    _structural_discipline->add_volume_load(0, *_flux_load);
    
    
    // element and material properties
    _parameters.resize(_n_vars);
    _parameter_functions.resize(_n_vars);
    
    // property values
    _E.reset(new MAST::Parameter("E", infile("youngs_modulus", 72.e9)));
    _nu.reset(new MAST::Parameter("nu", infile("poisson_ratio", 0.33)));
    _rho.reset(new MAST::Parameter("rho", infile("density", 2700.)));
    _kappa.reset(new MAST::Parameter("kappa", infile("shear_corr_factor", 5./6.)));
    
    // functions for property values
    _E_fn.reset(new MAST::ConstantFieldFunction("E", *_E));
    _nu_fn.reset(new MAST::ConstantFieldFunction("nu", *_nu));
    _rho_fn.reset(new MAST::ConstantFieldFunction("rho", *_rho));
    _kappa_fn.reset(new MAST::ConstantFieldFunction("kappa", *_kappa));
    
    
    _materials.reset(new MAST::IsotropicMaterialPropertyCard);
    
    const Real
    x0 = infile("x_div_loc", 0., 0),  // panel LE
    x1 = infile("x_div_loc", 0., 1);  // panel TE
    Real dx = (x1-x0)/(_n_stations-1);
    
    MAST::MaterialPropertyCardBase& mat = *_materials;
    // add the properties to the cards
    mat.add(*_E_fn);
    mat.add(*_nu_fn);
    mat.add(*_rho_fn);
    mat.add(*_kappa);
    
    // create the thickness variables
    for (unsigned int i=0; i<_n_stations; i++) {
        std::ostringstream oss;
        oss << "h_y_" << i;
        
        // now we need a parameter that defines the thickness at the
        // specified station and a constant function that defines the
        // field function at that location.
        MAST::Parameter* h_y               =
        new MAST::Parameter(oss.str(), infile("thickness", 0.002));
        
        MAST::ConstantFieldFunction* h_y_f =
        new MAST::ConstantFieldFunction("hy", *h_y);
        
        // add this to the density map
        _h_y_station_vals.insert(std::pair<Real, MAST::FieldFunction<Real>*>
                                 (x0+i*dx, h_y_f));
        
        // add the function to the parameter set
        _parameters[i]          = h_y->ptr();
        _parameter_functions[i] = h_y_f;
        
        // tell the assembly system about the sensitvity parameter
        _structural_discipline->add_parameter(*h_y);
    }
    
    // now create the h_y function and give it to the property card
    _h_y_fn.reset(new MAST::MultilinearInterpolation("hy", _h_y_station_vals));
    
    // create the variables for the h_z thickness
    _h_z.reset(new MAST::Parameter("hz", infile("width", 0.002)));
    _h_z_fn.reset(new MAST::ConstantFieldFunction("hz",*_h_z));
    _offset_h_y_fn.reset(new MAST::BeamOffset("hy_off", _h_y_fn->clone().release()));
    _offset_h_z.reset(new MAST::Parameter("hz_offset", 0.));
    _offset_h_z_fn.reset(new MAST::ConstantFieldFunction("hz_off", *_offset_h_z));
    
    MAST::Solid1DSectionElementPropertyCard *p =
    new MAST::Solid1DSectionElementPropertyCard;
    p->add(*_h_y_fn); // thickness
    p->add(*_h_z_fn); // width
    p->add(*_offset_h_y_fn); // thickness offset
    p->add(*_offset_h_z_fn); // width offset
    p->y_vector()(1) = 1.; // x-vector along x, y along y
    
    p->set_material(mat);
    p->set_strain(MAST::VON_KARMAN_STRAIN);
    _elem_properties.reset(p);
    p->init();
    
    _structural_discipline->set_property_for_subdomain(0, *_elem_properties);
    
    // create the function to calculate weight
    _weight = new MAST::Weight(*_structural_discipline);
    
    // create the mesh function to calculate the displacement
    std::vector<unsigned int> vars(1);
    vars[0] = _structural_sys->vars()[1]; // variable id for uy
    _disp_function.reset(new libMesh::MeshFunction(_static_system->get_equation_systems(),
                                                   *_static_system->solution,
                                                   _static_system->get_dof_map(),
                                                   vars));
    _disp_function->init();
    
    _disp_function_sens.resize(_n_vars);
    for (unsigned int i=0; i<_n_vars; i++) {
        _static_system->add_sensitivity_solution(i);
        _disp_function_sens[i] =
        new libMesh::MeshFunction(_static_system->get_equation_systems(),
                                  _static_system->get_sensitivity_solution(i),
                                  _static_system->get_dof_map(),
                                  vars);
        _disp_function_sens[i]->init();
    }
}




void
MAST::BeamBendingSizingOptimization::output(unsigned int iter,
                                            const std::vector<Real>& x,
                                            Real obj,
                                            const std::vector<Real>& fval,
                                            bool if_write_to_optim_file) const {
    
    libmesh_assert_equal_to(x.size(), _n_vars);
    
    // set the parameter values equal to the DV value
    for (unsigned int i=0; i<_n_vars; i++)
        *_parameters[i] = x[i];
    
    MAST::FunctionEvaluation::output(iter, x, obj, fval, if_write_to_optim_file);
}


#include <unistd.h>

int
main_4(int argc, char* const argv[]) {
    
    libMesh::LibMeshInit init(argc, argv);
    
    GetPot infile("input.in");
    
    std::ofstream output;
    output.open("optimization_output.txt", std::ofstream::out);
    
    MAST::GCMMAOptimizationInterface gcmma;
    
    // create and attach sizing optimization object
    MAST::BeamBendingSizingOptimization func_eval(init, output);
    
    // attach and optimize
    gcmma.attach_function_evaluation_object(func_eval);
    gcmma.optimize();
    
    output.close();
    
    return 0;
}


