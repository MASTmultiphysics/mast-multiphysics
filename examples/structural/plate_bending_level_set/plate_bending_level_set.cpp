/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2018  Manav Bhatia
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
#include "examples/structural/plate_bending_level_set/plate_bending_level_set.h"
#include "elasticity/structural_system_initialization.h"
#include "elasticity/structural_element_base.h"
#include "level_set/level_set_nonlinear_implicit_assembly.h"
#include "elasticity/structural_nonlinear_assembly.h"
#include "elasticity/structural_discipline.h"
#include "elasticity/stress_output_base.h"
#include "elasticity/structural_near_null_vector_space.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "property_cards/solid_2d_section_element_property_card.h"
#include "property_cards/isotropic_material_property_card.h"
#include "boundary_condition/dirichlet_boundary_condition.h"
#include "base/nonlinear_system.h"


// libMesh includes
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/parameter_vector.h"
#include "libmesh/diff_solver.h"
#include "libmesh/nonlinear_solver.h"


extern libMesh::LibMeshInit* __init;


class Phi:
public MAST::FieldFunction<Real> {
  
public:
    Phi(Real l1,
        Real l2,
        Real rmin,
        Real rmax):
    MAST::FieldFunction<Real>("Phi"),
    _l1(l1),
    _l2(l2),
    _rmin(rmin),
    _rmax(rmax) {
        
        libmesh_assert_less(rmax, _l1*.5);
        libmesh_assert_less(rmax, _l2*.5);
        libmesh_assert_less(rmin,   rmax);
    }
    virtual ~Phi() {}
    virtual void operator()(const libMesh::Point& p,
                            const Real t,
                            Real& v) const {
        
        libmesh_assert_less_equal(t, 1);
        v =
        pow(p(0)-_l1*.5, 2) +
        pow(p(1)-_l2*.5, 2) -
        pow((_rmin + t*(_rmax-_rmin)), 2);
    }
protected:
    Real
    _l1,
    _l2,
    _rmin,
    _rmax;
};


MAST::PlateBendingLevelSet::PlateBendingLevelSet():
_initialized(false) {

}

#include "plplot/plplot.h"
#include "level_set/level_set_intersection.h"
#include "level_set/sub_cell_fe.h"
#include <numeric>

void plot_elem(const libMesh::Elem& e) {
    
    unsigned int
    n  = e.n_nodes();
    
    RealVectorX
    x  = RealVectorX::Zero(n+1),
    y  = RealVectorX::Zero(n+1);
    
    for (unsigned int i=0; i<n+1; i++) {
        x(i) = e.point(i%n)(0);
        y(i) = e.point(i%n)(1);
    }
    
    plline(n+1, x.data(), y.data());
}


void plot_points(const std::vector<libMesh::Point>& pts) {

    unsigned int
    n  = pts.size();
    
    RealVectorX
    x  = RealVectorX::Zero(n),
    y  = RealVectorX::Zero(n);
    
    for (unsigned int i=0; i<n; i++) {
        x(i) = pts[i](0);
        y(i) = pts[i](1);
    }
    
    //char s[] = ".";
    //plstring(n, x.data(), y.data(), s);
    plpoin(n, x.data(), y.data(), -1);
}


void
MAST::PlateBendingLevelSet::init(libMesh::ElemType e_type,
                                 bool if_vk) {
    
    
    libmesh_assert(!_initialized);
    
    
    // length of domain
    _length     = 0.3,
    _width      = 0.3;
    
    
    // create the mesh
    _mesh       = new libMesh::ParallelMesh(__init->comm());
    
    // initialize the mesh with one element
    libMesh::MeshTools::Generation::build_square(*_mesh,
                                                 120, 120,
                                                 0, _length,
                                                 0, _width,
                                                 e_type);

    // create the equation system
    _eq_sys    = new  libMesh::EquationSystems(*_mesh);
    
    // create the libmesh system
    _sys       = &(_eq_sys->add_system<MAST::NonlinearSystem>("structural"));
    _phi_sys   = &(_eq_sys->add_system<libMesh::ExplicitSystem>("phi"));
    
    // FEType to initialize the system
    libMesh::FEType fetype (libMesh::FIRST, libMesh::LAGRANGE);
    
    // initialize the system to the right set of variables
    _structural_sys = new MAST::StructuralSystemInitialization(*_sys,
                                                               _sys->name(),
                                                               fetype);
    _discipline     = new MAST::StructuralDiscipline(*_eq_sys);
    
    _phi_sys->add_variable("phi", fetype);
    
    
//    {
//        Phi phi(_length, _width, _length*0.5*.1, _length*0.5*0.9);
//        libMesh::MeshBase::const_element_iterator
//        e_it   = _mesh->elements_begin(),
//        e_end  = _mesh->elements_end();
//
//        plsdev("xwin");
//        plinit();
//        plenv(0,.3,0,.3,0,0);
//
//        // plot the level set function
//        unsigned int npts = 100;
//        RealVectorX
//        phix = RealVectorX::Zero(npts+1),
//        phiy = RealVectorX::Zero(npts+1);
//        Real
//        pi = acos(-1.);
//
//        for (unsigned int i=0; i<npts+1; i++) {
//            phix(i) = 0.1*cos((1.*(i%npts))/(1.*npts)*2*pi)+0.15;
//            phiy(i) = 0.1*sin((1.*(i%npts))/(1.*npts)*2*pi)+0.15;
//        }
//        plline(npts+1, phix.data(), phiy.data());
//        plflush();
//
//        for ( ; e_it != e_end; e_it++) {
//
//
//            const libMesh::Elem *e = *e_it;
//            plcol0(1); // red color for elements
//            plot_elem(*e);
//            plflush();
//
//            {
//                MAST::FEBase fe(*_structural_sys);
//                fe.init(*e);
//                const std::vector<Real>& JxW = fe.get_JxW();
//                std::cout << "==== original JxW: " << std::accumulate(JxW.begin(),
//                                                                 JxW.end(), 0.) << std::endl;
//            }
//
//            MAST::LevelSetIntersection intersect;
//            intersect.init(phi, *e, 0.);
//
//            // now get elements on either side and plot
//            const std::vector<const libMesh::Elem *> &
//                    elems_low = intersect.get_sub_elems_negative_phi(),
//                    elems_hi = intersect.get_sub_elems_positive_phi();
//
//            plcol0(15); // white color for sub elements
//
//            for (unsigned int i = 0; i < elems_low.size(); i++) {
//                plot_elem(*elems_low[i]);
//                plflush();
//                // create FE
//                MAST::SubCellFE fe(*_structural_sys, intersect);
//                fe.init(*elems_low[i]);
//                const std::vector<Real>& JxW = fe.get_JxW();
//                const std::vector<libMesh::Point>& xyz = fe.get_xyz();
//                std::cout << "low: JxW: " << std::accumulate(JxW.begin(),
//                                                             JxW.end(), 0.) << std::endl;
//                plot_points(xyz);
//                plflush();
//
//                elems_low[i]->print_info();
//                for (unsigned int j=0; j<xyz.size(); j++)
//                    xyz[j].print();
//           }
//
//            for (unsigned int i=0; i<elems_hi.size(); i++) {
//                plot_elem(*elems_hi[i]);
//                plflush();
//
//                // create FE
//                MAST::SubCellFE fe(*_structural_sys, intersect);
//                fe.init(*elems_hi[i]);
//                const std::vector<Real>& JxW = fe.get_JxW();
//                const std::vector<libMesh::Point>& xyz = fe.get_xyz();
//                std::cout << "hi: JxW: " << std::accumulate(JxW.begin(),
//                                                             JxW.end(), 0.) << std::endl;
//                plot_points(fe.get_xyz());
//                plflush();
//
//                elems_hi[i]->print_info();
//                for (unsigned int j=0; j<xyz.size(); j++)
//                    xyz[j].print();
//            }
//
//        }
//        plend();
//    }
	
    // create and add the boundary condition and loads
    _dirichlet_bottom = new MAST::DirichletBoundaryCondition;
    _dirichlet_right  = new MAST::DirichletBoundaryCondition;
    _dirichlet_top    = new MAST::DirichletBoundaryCondition;
    _dirichlet_left   = new MAST::DirichletBoundaryCondition;
    
    _dirichlet_bottom->init (0, _structural_sys->vars());
    _dirichlet_right->init  (1, _structural_sys->vars());
    _dirichlet_top->init    (2, _structural_sys->vars());
    _dirichlet_left->init   (3, _structural_sys->vars());
    
    _discipline->add_dirichlet_bc(0, *_dirichlet_bottom);
    _discipline->add_dirichlet_bc(1,  *_dirichlet_right);
    _discipline->add_dirichlet_bc(2,    *_dirichlet_top);
    _discipline->add_dirichlet_bc(3,   *_dirichlet_left);
    _discipline->init_system_dirichlet_bc(*_sys);
    
    // initialize the equation system
    _eq_sys->init();
    
    // create the property functions and add them to the
    
    _th              = new MAST::Parameter("th",  0.01);
    _E               = new MAST::Parameter("E",  72.e9);
    _nu              = new MAST::Parameter("nu",   0.3);
    _kappa           = new MAST::Parameter("kappa",  5./6.);
    _zero            = new MAST::Parameter("zero",  0.);
    _press           = new MAST::Parameter( "p",  3.e7);
    
    
    
    // prepare the vector of parameters with respect to which the sensitivity
    // needs to be benchmarked
    _params_for_sensitivity.push_back(_E);
    _params_for_sensitivity.push_back(_nu);
    _params_for_sensitivity.push_back(_th);
    
    
    
    _th_f            = new MAST::ConstantFieldFunction("h",           *_th);
    _E_f             = new MAST::ConstantFieldFunction("E",            *_E);
    _nu_f            = new MAST::ConstantFieldFunction("nu",          *_nu);
    _kappa_f         = new MAST::ConstantFieldFunction("kappa",    *_kappa);
    _hoff_f          = new MAST::ConstantFieldFunction("off",       *_zero);
    _press_f         = new MAST::ConstantFieldFunction("pressure", *_press);
    
    // initialize the load
    _p_load          = new MAST::BoundaryConditionBase(MAST::SURFACE_PRESSURE);
    _p_load->add(*_press_f);
    _discipline->add_volume_load(0, *_p_load);
    
    // create the material property card
    _m_card         = new MAST::IsotropicMaterialPropertyCard;
    
    // add the material properties to the card
    _m_card->add(*_E_f);
    _m_card->add(*_nu_f);
    _m_card->add(*_kappa_f);
    
    // create the element property card
    _p_card         = new MAST::Solid2DSectionElementPropertyCard;
    
    // add the section properties to the card
    _p_card->add(*_th_f);
    _p_card->add(*_hoff_f);
    
    // tell the section property about the material property
    _p_card->set_material(*_m_card);
    _p_card->set_bending_model(MAST::MINDLIN);
    if (if_vk) _p_card->set_strain(MAST::VON_KARMAN_STRAIN);
    
    _discipline->set_property_for_subdomain(0, *_p_card);
    
    
    // create the output objects, one for each element
    libMesh::MeshBase::const_element_iterator
    e_it    = _mesh->elements_begin(),
    e_end   = _mesh->elements_end();
    
    // points where stress is evaluated
    std::vector<libMesh::Point> pts;
    if (e_type == libMesh::QUAD4 ||
        e_type == libMesh::QUAD8 ||
        e_type == libMesh::QUAD9) {
        
        pts.push_back(libMesh::Point(-1/sqrt(3), -1/sqrt(3), 1.)); // upper skin
        pts.push_back(libMesh::Point(-1/sqrt(3), -1/sqrt(3),-1.)); // lower skin
        pts.push_back(libMesh::Point( 1/sqrt(3), -1/sqrt(3), 1.)); // upper skin
        pts.push_back(libMesh::Point( 1/sqrt(3), -1/sqrt(3),-1.)); // lower skin
        pts.push_back(libMesh::Point( 1/sqrt(3),  1/sqrt(3), 1.)); // upper skin
        pts.push_back(libMesh::Point( 1/sqrt(3),  1/sqrt(3),-1.)); // lower skin
        pts.push_back(libMesh::Point(-1/sqrt(3),  1/sqrt(3), 1.)); // upper skin
        pts.push_back(libMesh::Point(-1/sqrt(3),  1/sqrt(3),-1.)); // lower skin
    }
    else if (e_type == libMesh::TRI3 ||
             e_type == libMesh::TRI6) {
        
        pts.push_back(libMesh::Point(1./3., 1./3., 1.)); // upper skin
        pts.push_back(libMesh::Point(1./3., 1./3.,-1.)); // lower skin
        pts.push_back(libMesh::Point(2./3., 1./3., 1.)); // upper skin
        pts.push_back(libMesh::Point(2./3., 1./3.,-1.)); // lower skin
        pts.push_back(libMesh::Point(1./3., 2./3., 1.)); // upper skin
        pts.push_back(libMesh::Point(1./3., 2./3.,-1.)); // lower skin
    }
    else
        libmesh_assert(false); // should not get here
    
    for ( ; e_it != e_end; e_it++) {
        
        MAST::StressStrainOutputBase * output = new MAST::StressStrainOutputBase;
        
        // tell the object to evaluate the data for this object only
        std::set<const libMesh::Elem*> e_set;
        e_set.insert(*e_it);
        output->set_elements_in_domain(e_set);
        output->set_points_for_evaluation(pts);
        output->set_volume_loads(_discipline->volume_loads());
        _outputs.push_back(output);
        
        _discipline->add_volume_output((*e_it)->subdomain_id(), *output);
    }
    
    _initialized = true;
}







MAST::PlateBendingLevelSet::~PlateBendingLevelSet() {
    
    if (_initialized) {
        
        delete _m_card;
        delete _p_card;
        
        delete _p_load;
        delete _dirichlet_bottom;
        delete _dirichlet_right;
        delete _dirichlet_top;
        delete _dirichlet_left;
        
        delete _th_f;
        delete _E_f;
        delete _nu_f;
        delete _kappa_f;
        delete _hoff_f;
        delete _press_f;
        
        delete _th;
        delete _E;
        delete _nu;
        delete _kappa;
        delete _zero;
        delete _press;
        
        
        
        delete _eq_sys;
        delete _mesh;
        
        delete _discipline;
        delete _structural_sys;
        
        // iterate over the output quantities and delete them
        std::vector<MAST::StressStrainOutputBase*>::iterator
        it   =   _outputs.begin(),
        end  =   _outputs.end();
        
        for ( ; it != end; it++)
            delete *it;
        
        _outputs.clear();
    }
}



MAST::Parameter*
MAST::PlateBendingLevelSet::get_parameter(const std::string &nm) {
    
    libmesh_assert(_initialized);
    
    MAST::Parameter *rval = nullptr;
    
    // look through the vector of parameters to see if the name is available
    std::vector<MAST::Parameter*>::iterator
    it   =  _params_for_sensitivity.begin(),
    end  =  _params_for_sensitivity.end();
    
    bool
    found = false;
    
    for ( ; it != end; it++) {
        
        if (nm == (*it)->name()) {
            rval    = *it;
            found   = true;
        }
    }
    
    // if the param was not found, then print the message
    if (!found) {
        libMesh::out
        << std::endl
        << "Parameter not found by name: " << nm << std::endl
        << "Valid names are: "
        << std::endl;
        for (it = _params_for_sensitivity.begin(); it != end; it++)
            libMesh::out << "   " << (*it)->name() << std::endl;
        libMesh::out << std::endl;
    }
    
    return rval;
}


void
_project_phi(MAST::FieldFunction<Real>& phi,
             libMesh::ExplicitSystem& sys) {
    
    // now create a function and use it for initialization
    class SolutionFunction:
    public libMesh::FunctionBase<Real> {
    public:
        SolutionFunction(MAST::FieldFunction<Real>& phi):
        libMesh::FunctionBase<Real>(),
        _phi(phi) { }
        
        virtual std::unique_ptr<libMesh::FunctionBase<Real> > clone () const {
            libMesh::FunctionBase<Real> *rval = new SolutionFunction(_phi);
            return std::unique_ptr<libMesh::FunctionBase<Real> >(rval);
        }
        
        // this should not get called
        virtual Real operator()
        (const libMesh::Point& p, const Real time) {libmesh_assert(false);}
        
        virtual void
        operator() (const libMesh::Point& p,
                    const Real time,
                    libMesh::DenseVector<Real>& output) {
            Real v;
            _phi(p, time, v);
            output(0) = v;
        }
    protected:
        MAST::FieldFunction<Real>& _phi;
    };
    
    SolutionFunction sol_func(phi);
    
    sys.project_solution(&sol_func);
}



const libMesh::NumericVector<Real>&
MAST::PlateBendingLevelSet::solve(bool if_write_output) {
    
    
    libmesh_assert(_initialized);
    
    bool if_vk = (_p_card->strain_type() == MAST::VON_KARMAN_STRAIN);
    
    // set the number of load steps
    unsigned int
    n_steps = 10;
    //if (if_vk) n_steps = 50;
    
    Real
    p0      = (*_press)();
    
    Phi phi(_length, _width, _length*0.5*.2, _length*0.5*0.9);
    
    // create the nonlinear assembly object
    MAST::LevelSetNonlinearImplicitAssembly           assembly;
    MAST::StructuralNonlinearAssemblyElemOperations   elem_ops;
    
    assembly.attach_discipline_and_system(elem_ops,
                                          *_discipline,
                                          *_structural_sys,
                                          phi);
    
    MAST::NonlinearSystem& nonlin_sys = assembly.system();
    
    // zero the solution before solving
    nonlin_sys.solution->zero();
    this->clear_stresss();
    
    MAST::StructuralNearNullVectorSpace nsp;
    nonlin_sys.nonlinear_solver->nearnullspace_object = &nsp;
    
    
    libMesh::ExodusII_IO exodus_writer(*_mesh);
    
    // now iterate over the load steps
    for (unsigned int i=0; i<n_steps; i++) {
        
        //(*_press)()  =  p0*(i+1.)/(1.*n_steps);
        nonlin_sys.time = (i+1.)/(1.*n_steps);
        _phi_sys->time  = (i+1.)/(1.*n_steps);
        
        // reinit constraints for the new level set
        nonlin_sys.reinit_constraints();
        _project_phi(phi, *_phi_sys);
        
        libMesh::out
        << "Load step: " << i << "  : t = " << nonlin_sys.time << std::endl;
        
        nonlin_sys.solve();
        
        // evaluate the outputs
        this->clear_stresss();
        assembly.calculate_outputs(*(_sys->solution));
        
        
        if (if_write_output) {
            
            libMesh::out << "Writing output to : output.exo" << std::endl;
            
            // write the solution for visualization
            _discipline->update_stress_strain_data();
            exodus_writer.write_timestep("output.exo",
                                         *_eq_sys,
                                         i+1,
                                         (1.*i)/(1.*(n_steps-1)));
        }
    }
    assembly.clear_discipline_and_system();
    
    return *(_sys->solution);
}





const libMesh::NumericVector<Real>&
MAST::PlateBendingLevelSet::sensitivity_solve(MAST::Parameter& p,
                                              bool if_write_output) {
    
    libmesh_assert(_initialized);
    
    _discipline->add_parameter(p);

    Phi phi(_length, _width, _length*0.5*.1, _length*0.5*0.9);

    // create the nonlinear assembly object
    MAST::LevelSetNonlinearImplicitAssembly           assembly;
    MAST::StructuralNonlinearAssemblyElemOperations   elem_ops;
    
    assembly.attach_discipline_and_system(elem_ops,
                                          *_discipline,
                                          *_structural_sys,
                                          phi);

    MAST::NonlinearSystem& nonlin_sys = assembly.system();
    
    libMesh::ParameterVector params;
    params.resize(1);
    params[0]  =  p.ptr();
    
    // zero the solution before solving
    nonlin_sys.add_sensitivity_solution(0).zero();
    this->clear_stresss();
    
    nonlin_sys.sensitivity_solve(params);
    
    // evaluate sensitivity of the outputs
    assembly.calculate_output_sensitivity(params,
                                          true,    // true for total sensitivity
                                          *(_sys->solution));
    
    
    assembly.clear_discipline_and_system();
    _discipline->remove_parameter(p);
    
    // write the solution for visualization
    if (if_write_output) {
        
        std::ostringstream oss1, oss2;
        oss1 << "output_"        << p.name() << ".exo";
        oss2 << "stress_output_" << p.name() << ".exo";
        
        libMesh::out
        << "Writing sensitivity output to : " << oss1.str()
        << "  and stress/strain sensitivity to : " << oss2.str()
        << std::endl;
        
        
        _sys->solution->swap(_sys->get_sensitivity_solution(0));
        
        // write the solution for visualization
        _discipline->update_stress_strain_data( &p);
        libMesh::ExodusII_IO(*_mesh).write_equation_systems(oss1.str(),
                                                            *_eq_sys);
        
        _sys->solution->swap(_sys->get_sensitivity_solution(0));
    }
    
    return _sys->get_sensitivity_solution(0);
}



void
MAST::PlateBendingLevelSet::clear_stresss() {
    
    libmesh_assert(_initialized);
    
    // iterate over the output quantities and delete them
    std::vector<MAST::StressStrainOutputBase*>::iterator
    it   =   _outputs.begin(),
    end  =   _outputs.end();
    
    for ( ; it != end; it++)
        (*it)->clear(false);
}



