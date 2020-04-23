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

// C/C++ includes.
#include <iostream>

// MAST includes.
#include "examples/fluid/meshing/cylinder.h"
#include "examples/fluid/meshing/naca0012.h"
#include "examples/fluid/meshing/panel_mesh_2D.h"
#include "examples/fluid/meshing/panel_mesh_3D.h"
#include "examples/fluid/meshing/naca0012_wing.h"
#include "examples/base/input_wrapper.h"
#include "base/nonlinear_system.h"
#include "base/transient_assembly.h"
#include "base/boundary_condition_base.h"
#include "base/field_function_base.h"
#include "base/parameter.h"
#include "boundary_condition/dirichlet_boundary_condition.h"
#include "fluid/conservative_fluid_system_initialization.h"
#include "fluid/conservative_fluid_discipline.h"
#include "fluid/conservative_fluid_transient_assembly.h"
#include "fluid/flight_condition.h"
#include "fluid/integrated_force_output.h"
#include "solver/first_order_newmark_transient_solver.h"
#include "solver/stabilized_first_order_transient_sensitivity_solver.h"

// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_elem_type.h"    // ElemType
#include "libmesh/fe_type.h"           // FEFamily, Order
#include "libmesh/parallel_mesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/dof_map.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/nonlinear_solver.h"
#include "libmesh/periodic_boundary.h"

/*!
 *    class computes the far-field solutions in zones 1, 2, 3 based on the
 *    initial shock angle \p beta. Note, that \p beta is the same as
 *    \p sigma in Nathan, et al's 2018 paper. The relations are based on
 *    standard theory concerning oblique shocks.
 */
class ShockImpingementSolutions:
public MAST::FieldFunction<RealVectorX> {
public:
    
    /*!
     *    Initializes the far-field solutions in zones 1, 2, 3 based on the
     *    initial shock angle \p beta. Note, that \p beta is the same as
     *    \p sigma in Nathan, et al's 2018 paper. \p cp, \p cv are the
     *    coefficients at constant pressure and volume. \p M1 is the Mach
     *    number in zone 1 (before shock), \p rho1 is the fluid density in
     *    zone 1, \p T1 is the fluid temprature (absolute value) in zone 1,
     *    \p beta is the shock angle. \p x is the physical location along
     *    the x-axis where the shock impinges the panel.
     */
    ShockImpingementSolutions(Real cp,
                              Real cv,
                              Real M1,
                              Real rho1,
                              Real T1,
                              Real beta,
                              Real x0):
    MAST::FieldFunction<RealVectorX>("fluid_solution"),
    _x0      (x0),
    _beta1   (beta),
    _beta2   (0.),
    _theta1  (0.) {
        /*gamma     = 1.4;
         R         = 1003-716;
         M1        = 2.0;
         rho1      = 1.35;
         T1        = 300;
         beta_deg  = 34.799;
         beta      = beta_deg * pi / 180.;
         */
        
        libmesh_assert_greater(M1, 1.);
        
        Real
        gamma       = cp/cv,
        R           = cp-cv,
        p1          = rho1 * R * T1,
        M2          = _M2(gamma, beta, M1),
        rho2        = _rho2(rho1, gamma, beta, M1),
        p2          = _p2(p1, gamma, beta, M1),
        T2          = p2/rho2/R;
        
        _theta1     = _theta(gamma, beta, M1);
        _beta2      = _find_beta(gamma, _theta1, beta, M2);
        
        // expected theta after shock is the same as first theta. Now we find beta that gives us the theta for that
        Real
        M3          = _M2(gamma, _beta2, M2),
        rho3        = _rho2(rho2, gamma, _beta2, M2),
        p3          = _p2(p2, gamma, _beta2, M2),
        T3          = p3/rho3/R,
        a           = 0.,
        u1          = 0.,
        u2          = 0.;
        libMesh::out << "p3/p1: " << p3/p1 << std::endl;
        
        _sol1 = RealVectorX::Zero(4);
        _sol2 = RealVectorX::Zero(4);
        _sol3 = RealVectorX::Zero(4);
        
        // zone 1
        a        = sqrt(gamma*R*T1);
        u1       = a*M1;
        u2       = 0.;
        _sol1(0)  = rho1;
        _sol1(1)  = rho1 * u1;
        _sol1(2)  = rho1 * u2;
        _sol1(3)  = rho1 * (cv*T1 + .5*(u1*u1 + u2*u2));
        
        // zone 2
        a        =  sqrt(gamma*R*T2);
        u1       =  a*M2*cos(_theta1);
        u2       = -a*M2*sin(_theta1);
        _sol2(0)  = rho2;
        _sol2(1)  = rho2 * u1;
        _sol2(2)  = rho2 * u2;
        _sol2(3)  = rho2 * (cv*T2 + .5*(u1*u1 + u2*u2));
        
        // zone 3
        a        = sqrt(gamma*R*T3);
        u1       = a*M3;
        u2       = 0.;
        _sol3(0)  = rho3;
        _sol3(1)  = rho3 * u1;
        _sol3(2)  = rho3 * u2;
        _sol3(3)  = rho3 * (cv*T3 + .5*(u1*u1 + u2*u2));
    }
    
    /*!
     *    identifies the far-field solution based on spatial location \p p
     */
    virtual void operator() (const libMesh::Point& p,
                             const Real t,
                             RealVectorX& v) const {
        
        if (p(0) <= _x0) {
            if (p(1) <= tan(_beta1)*(_x0-p(0)))
                v = _sol1;
            else
                v = _sol2;
        }
        else {
            if (p(1) <= tan(_beta2-_theta1)*(p(0)-_x0))
                v = _sol3;
            else
                v = _sol2;
        }
    }
    
protected:
    
    Real _theta(Real gamma, Real beta, Real M1) {
        return atan( 2./tan(beta) * ((pow(M1*sin(beta),2) - 1.) / (pow(M1,2)*(gamma+cos(2.*beta))+2)) );
    }
    
    Real _p2(Real p1, Real gamma, Real beta, Real M1) {
        return p1 * (1. + 2*gamma/(gamma+1.) * (pow(M1*sin(beta),2) -1));
    }
    
    Real _rho2(Real rho1, Real gamma, Real beta, Real M1) {
        return rho1 * ((gamma+1.)*pow(M1*sin(beta),2)/ ((gamma-1.) * pow(M1*sin(beta),2) + 2.));
    }
    
    Real _M2(Real gamma, Real beta, Real M1) {
        
        Real thetav = _theta(gamma, beta, M1);
        return 1./(sin(beta-thetav)) *
        sqrt((1.+(gamma-1.)/2. * pow(M1*sin(beta),2))/(gamma*pow(M1*sin(beta),2)-(gamma-1.)/2.));
    }
    
    Real _find_beta(Real gamma, Real theta0, Real betav0, Real M1) {
        
        bool cont = true;
        Real
        betav    = betav0,
        delta    = 1.e-4,
        tol      = 1.e-6,
        res      = 0.,
        dbeta    = 0.,
        jac      = 0.;
        
        while (cont) {
            res        = _theta(gamma, betav, M1) - theta0;
            if (fabs(res) > tol) {
                
                jac     = (_theta(gamma, betav+delta, M1)- _theta(gamma, betav, M1))/delta;
                dbeta   = -res/jac;
                betav   = betav + dbeta;
            }
            else
                cont = false;
        }
        
        return betav;
    }
    
    Real _x0, _beta1, _beta2, _theta1;
    RealVectorX _sol1, _sol2, _sol3;
};


//
// BEGIN_TRANSLATE Flow analysis
//
//   \tableofcontents
//
//
//   This example computes the steady-state flow about rigid shapes.
//   A Newmark time integration is used for time-stepping. Since a time-accurate
//   flow solution is not sought the Newton-Raphson scheme at each time-step
//   is only partially converged. The following options can be used to
//   tune the SNES solver in PETSc to accomplish this:
//   ```
//   -snes_linesearch_type basic -snes_linesearch_damping 1.0 -snes_linesearch_max_it 1 -snes_max_it 2
//   ```
//
//  This class stores all the data structures necessary to setup a flow
//  analysis
class FlowAnalysis {
    
protected:
    
    libMesh::LibMeshInit&                          _init;
    MAST::Examples::GetPotWrapper&                 _input;
    unsigned int                                   _dim;
    libMesh::UnstructuredMesh*                     _mesh;
    libMesh::EquationSystems*                      _eq_sys;
    MAST::NonlinearSystem*                         _sys;
    MAST::ConservativeFluidSystemInitialization*   _sys_init;
    MAST::ConservativeFluidDiscipline*             _discipline;
    MAST::FlightCondition*                         _flight_cond;
    
    libMesh::ExodusII_IO*                          _output;
    
    libMesh::FEType                                _fetype;

    std::set<MAST::BoundaryConditionBase*>         _boundary_conditions;
    std::set<libMesh::PeriodicBoundary*>           _periodic_boundaries;

    //
    // \section mesh Mesh generation
    // This will initialize the mesh
    void _init_mesh(bool mesh, bool bc) {
        
        // The mesh is created using classes written in MAST. The particular
        // mesh to be used can be selected using the input parameter
        // ` mesh=val `, where `val` can be one of the following:
        //   - `naca0012` for flow over a NACA 0012 airfoil
        //   - `cylinder` for flow over a cylinder
        //   - `naca0012_wing` for flow over swept wing with NACA0012 section
        //   - `panel_2D` for 2D flow analysis over a panel
        //   - `panel_3D` for 3D flow analysis over a panel
        //   - `mfu`      for 3D flow analysis over a minimal flow unit
        //
        // The meshing and boundary conditions for each flow analysis case
        // can be modified using suitable parameters, which are described in this
        // example.
        std::string
        s  = _input("mesh",
                    "type of mesh to be analyzed {naca0012, cylinder, naca0012_wing, panel_2D, panel_3D, mfu}",
                    "naca0012");

        if (s == "naca0012")
            _init_naca0012(mesh, bc);
        else if (s == "cylinder")
            _init_cylinder(mesh, bc);
        else if (s == "naca0012_wing")
            _init_naca0012_wing(mesh, bc);
        else if (s == "panel_2d")
            _init_panel_2D(mesh, bc);
        else if (s == "panel_3d")
            _init_panel_3D(mesh, bc);
        else if (s == "mfu")
            _init_minimal_flow_unit(mesh, bc);
        else
            libmesh_error_msg("unknown mesh type");
    }
    
    // \subsection naca0012_mesh NACA0012
    // If \p mesh is \p true then the mesh will be initialized. If \p bc
    // is true then the boundary conditions will be initialized
    void _init_naca0012(bool mesh, bool bc) {
        
        if (mesh) {
            
            _dim                 = 2;

            const unsigned int
            radial_divs         = _input("n_radial_elems", "number of elements in the radial direction from airfoil to far-field", 20),
            quarter_divs        = _input("n_quarter_elems", "number of elements in the quarter arc along the circumferencial direction", 20);
            
            const Real
            r                   = _input("chord", "chord of the airfoil", 0.1),
            l_by_r              = _input("l_by_r", "far-field distance to airfoil radius ratio", 5.),
            h_ff_by_h_r         = _input("h_far_field_by_h_r", "relative element size at far-field boundary to element size at airfoil", 50.);
            
            std::string
            t = _input("elem_type", "type of geometric element in the mesh", "quad4");
            
            libMesh::ElemType
            e_type = libMesh::Utility::string_to_enum<libMesh::ElemType>(t);
            
            // initialize the mesh
            MAST::Examples::NACA0012Mesh2D().mesh(r,
                                                  r*l_by_r,
                                                  radial_divs,
                                                  quarter_divs,
                                                  h_ff_by_h_r,
                                                  *_mesh,
                                                  e_type);
            
        }
        
        if (bc) {
            
            libmesh_assert(_flight_cond);
            
            bool
            if_viscous = _flight_cond->gas_property.if_viscous;

            if (if_viscous) {

                MAST::DirichletBoundaryCondition
                *dirichlet = new MAST::DirichletBoundaryCondition;

                std::vector<unsigned int> constrained_vars(2);
                constrained_vars[0] = _sys_init->vars()[1];
                constrained_vars[1] = _sys_init->vars()[2];
                dirichlet->init (3, constrained_vars);
                _discipline->add_dirichlet_bc(3, *dirichlet);
                _discipline->init_system_dirichlet_bc(*_sys);

                _boundary_conditions.insert(dirichlet);
            }

            // create the boundary conditions for slip-wall, symmetry and far-field
            MAST::BoundaryConditionBase
            *far_field   = new MAST::BoundaryConditionBase(MAST::FAR_FIELD),
            // if a viscous analysis is requested then set the wall to be a no-slip
            // wall. Otherwise, use a slip wall for inviscid analysis
            *wall        = new MAST::BoundaryConditionBase(if_viscous?
                                                           MAST::NO_SLIP_WALL:
                                                           MAST::SLIP_WALL);
            
            // airfoil surface is boundary id 3 ...
            _discipline->add_side_load(   3, *wall);
            // boundary id 1 is far field
            _discipline->add_side_load(   1, *far_field);
            
            // store the pointers for later deletion in the destructor
            _boundary_conditions.insert(far_field);
            _boundary_conditions.insert(wall);
        }
    }

    // \subsection cylinder_mesh Cylinder
    void _init_cylinder(bool mesh, bool bc) {
        
        if (mesh) {

            _dim                 = 2;

            const unsigned int
            radial_divs         = _input("n_radial_elems", "number of elements in the radial direction from cylinder to far-field", 20),
            quarter_divs        = _input("n_quarter_elems", "number of elements in the quarter arc along the circumferencial direction", 20);
            
            const Real
            r                   = _input("radius", "radius of the cylinder", 0.1),
            l_by_r              = _input("l_by_r", "far-field distance to cylinder radius ratio", 5.),
            h_ff_by_h_r         = _input("h_far_field_by_h_r", "relative element size at far-field boundary to element size at cylinder", 50.);
            
            std::string
            t = _input("elem_type", "type of geometric element in the mesh", "quad4");
            
            libMesh::ElemType
            e_type = libMesh::Utility::string_to_enum<libMesh::ElemType>(t);
            
            // initialize the mesh
            MAST::Examples::CylinderMesh2D().mesh(r,
                                                  r*l_by_r,
                                                  radial_divs,
                                                  quarter_divs,
                                                  h_ff_by_h_r,
                                                  *_mesh,
                                                  e_type,
                                                  false, 0, 0);

        }
        
        if (bc) {

            libmesh_assert(_flight_cond);
            
            bool
            if_viscous = _flight_cond->gas_property.if_viscous;
            
            if (if_viscous) {
                
                MAST::DirichletBoundaryCondition
                *dirichlet = new MAST::DirichletBoundaryCondition;
                
                std::vector<unsigned int> constrained_vars(2);
                constrained_vars[0] = _sys_init->vars()[1];
                constrained_vars[1] = _sys_init->vars()[2];
                dirichlet->init (3, constrained_vars);
                _discipline->add_dirichlet_bc(3, *dirichlet);
                _discipline->init_system_dirichlet_bc(*_sys);
                
                _boundary_conditions.insert(dirichlet);
            }
            
            // create the boundary conditions for slip-wall, symmetry and far-field
            MAST::BoundaryConditionBase
            *far_field   = new MAST::BoundaryConditionBase(MAST::FAR_FIELD),
            // if a viscous analysis is requested then set the wall to be a no-slip
            // wall. Otherwise, use a slip wall for inviscid analysis
            *wall        = new MAST::BoundaryConditionBase(if_viscous?
                                                           MAST::NO_SLIP_WALL:
                                                           MAST::SLIP_WALL);
            
            // Cylinder surface is boundary id 3 ...
            _discipline->add_side_load(   3, *wall);
            // boundary id 1 is far field
            _discipline->add_side_load(   1, *far_field);

            // store the pointers for later deletion in the destructor
            _boundary_conditions.insert(far_field);
            _boundary_conditions.insert(wall);

        }
    }

    // \subsection naca0012_wing_mesh NACA0012 Wing
    void _init_naca0012_wing(bool mesh, bool bc) {
     
        if (mesh) {
            
            _dim  = 3;
            
            const unsigned int
            radial_divs_chord  = _input("radial_divs_chord", "number of elements along the radial direction from mid-chord to trailing-edge", 4),
            radial_divs_chord_to_farfield = _input("radial_divs_chord_to_farfield", "number of elements along the radial direction from trailing-edge to far-field", 20),
            quarter_divs        = _input("n_quarter_elems", "number of elements in the quarter arc along the circumferencial direction", 20),
            spanwise_divs       = _input("spanwise_divs", "number of elements along the span", 20),
            span_to_farfield_divs= _input("span_to_farfield_divs", "number of elements along the spanwise direction from wing-tip to far-field", 20);
            
            const Real
            root_chord          = _input("root_chord", "chord length at thw wing root", 0.5),
            taper_ratio         = _input("taper_ratio", "ratio of chord at wing-tip to chord at root", 0.5),
            mid_chord_sweep     = _input("mid_chord_sweep", "sweep of the mid-chord in radians", 0.35),
            far_field_radius_to_root_chord = _input("far_field_radius_to_root_chord", "radial-far field distance in multiples of wing chord", 15.),
            span                = _input("span", "wing span measured along y-axis", 2.),
            spanwise_farfield   = _input("spanwise_farfield", "distance of far-field boundary along y-axis from wing root", 10),
            radial_elem_size_ratio = _input("radial_elem_size_ratio", "ratio of element radial size at far-field to the size at wing surface", 10.),
            spanwise_elem_size_ratio = _input("spanwise_elem_size_ratio", "ratio of element spanwise size at far-field to the size at wing surface", 10.);
            
            std::string
            t = _input("elem_type", "type of geometric element in the mesh", "hex8");
            
            libMesh::ElemType
            e_type = libMesh::Utility::string_to_enum<libMesh::ElemType>(t);
            
            // initialize the mesh
            MAST::Examples::NACA0012WingMesh3D().mesh(root_chord,
                                                      taper_ratio,
                                                      mid_chord_sweep,
                                                      far_field_radius_to_root_chord,
                                                      span,
                                                      spanwise_farfield,
                                                      radial_divs_chord,
                                                      radial_divs_chord_to_farfield,
                                                      quarter_divs,
                                                      spanwise_divs,
                                                      span_to_farfield_divs,
                                                      radial_elem_size_ratio,
                                                      spanwise_elem_size_ratio,
                                                      *_mesh,
                                                      e_type);
            
        }
        
        if (bc) {
            
            libmesh_assert(_flight_cond);
            
            bool
            if_viscous = _flight_cond->gas_property.if_viscous;
            
            // create the boundary conditions for slip-wall, symmetry and far-field
            MAST::BoundaryConditionBase
            *far_field   = new MAST::BoundaryConditionBase(MAST::FAR_FIELD),
            *symmetry    = new MAST::BoundaryConditionBase(MAST::SYMMETRY_WALL),
            // if a viscous analysis is requested then set the wall to be a no-slip
            // wall. Otherwise, use a slip wall for inviscid analysis
            *wall        = new MAST::BoundaryConditionBase(if_viscous?
                                                           MAST::NO_SLIP_WALL:
                                                           MAST::SLIP_WALL);
            
            // wing surface is boundary id 4 ...
            _discipline->add_side_load(   4, *wall);      // wing surface
            _discipline->add_side_load(   0, *symmetry);  // root
            _discipline->add_side_load(   2, *far_field); // radial  far-field
            _discipline->add_side_load(   5, *far_field); // spanwise farfield

            // store the pointers for later deletion in the destructor
            _boundary_conditions.insert(far_field);
            _boundary_conditions.insert(symmetry);
            _boundary_conditions.insert(wall);
        }

    }

    // \subsection mfu_mesh Minimal Flow Unit
    void _init_minimal_flow_unit(bool mesh, bool bc) {
        
        if (mesh) {
            
            libmesh_assert(_sys);
            
            _dim  = 3;
            
            const unsigned int
            streamline_divs  = _input("streamline_divs", "number of elements along the streamwise direction", 20),
            crossflow_divs   = _input("crossflow_divs", "number of elements along the crossflow direction",   20),
            height_divs      = _input("height_divs", "number of elements along the height", 20);
            
            const Real
            pi               = acos(-1.),
            length           = _input("length", "streamwise length of flow unit", pi),
            width            = _input("width", "crossflow width of the flow unit", 0.34*pi),
            height           = _input("height", "distance between upper and lower walls", 2.);
            

            // initialize the periodic boundary conditions
            // libMesh numbering:
            //  0  -> back    -- along z
            //  1  -> bottom  -- along y
            //  2  -> right   -- along x
            //  3  -> top     -- along y
            //  4  -> left    -- along x
            //  5  -> front   -- along z
            
            libMesh::PeriodicBoundary
            *periodic_streamwise = new libMesh::PeriodicBoundary(libMesh::Point(length, 0., 0.)),
            *periodic_crossflow  = new libMesh::PeriodicBoundary(libMesh::Point(0.,  width, 0.));
            
            // left,right along x-axis
            periodic_streamwise->myboundary     = 4;
            periodic_streamwise->pairedboundary = 2;

            // bottom-top along y-axis
            periodic_crossflow->myboundary      = 1;
            periodic_crossflow->pairedboundary  = 3;

            _sys->get_dof_map().add_periodic_boundary(*periodic_streamwise);
            _sys->get_dof_map().add_periodic_boundary(*periodic_crossflow);

            
            
            std::string
            t = _input("elem_type", "type of geometric element in the mesh", "hex8");
            
            libMesh::ElemType
            e_type = libMesh::Utility::string_to_enum<libMesh::ElemType>(t);
            
            // initialize the mesh
            libMesh::MeshTools::Generation::build_cube(*_mesh,
                                                       streamline_divs,
                                                       crossflow_divs,
                                                       height_divs,
                                                       -length/2., length/2.,
                                                       -width/2.,   width/2.,
                                                       -height/2., height/2.,
                                                       e_type);
            
        }
        
        if (bc) {
            
            libmesh_assert(_flight_cond);
            
            // this assumes a viscous analysis
            libmesh_assert(_flight_cond->gas_property.if_viscous);
            
            // libMesh numbering:
            //  0  -> back    -- along z
            //  1  -> bottom  -- along y
            //  2  -> right   -- along x
            //  3  -> top     -- along y
            //  4  -> left    -- along x
            //  5  -> front   -- along z

            MAST::DirichletBoundaryCondition
            *dirichlet_bottom = new MAST::DirichletBoundaryCondition,
            *dirichlet_top    = new MAST::DirichletBoundaryCondition;
            
            std::vector<unsigned int> constrained_vars(3);
            constrained_vars[0] = _sys_init->vars()[1];
            constrained_vars[1] = _sys_init->vars()[2];
            constrained_vars[2] = _sys_init->vars()[3];
            dirichlet_bottom->init (0, constrained_vars);
            dirichlet_top->init    (5, constrained_vars);
            _discipline->add_dirichlet_bc(0, *dirichlet_bottom);
            _discipline->add_dirichlet_bc(5, *dirichlet_top);
            _discipline->init_system_dirichlet_bc(*_sys);
            
            // create the boundary conditions for slip-wall, symmetry and far-field
            MAST::BoundaryConditionBase
            *wall        = new MAST::BoundaryConditionBase(MAST::NO_SLIP_WALL);

            // libMesh defines back/front along the z-axis
            _discipline->add_side_load(   0, *wall);      // back
            _discipline->add_side_load(   5, *wall);      // front
            
            // store the pointers for later deletion in the destructor
            _boundary_conditions.insert(wall);
            _boundary_conditions.insert(dirichlet_top);
            _boundary_conditions.insert(dirichlet_bottom);
        }
        
    }

    // \subsection panel_2d_mesh Two-dimensional Panel Flow
    void _init_panel_2D(bool mesh, bool bc) {
        
        if (mesh) {
            
            _dim                 = 2;
            
            const unsigned int
            nx_divs             = 3,
            ny_divs             = 1,
            n_divs_ff_to_panel  = _input("n_divs_farfield_to_panel", "number of element divisions from far-field to panel", 30),
            n_divs_panel        = _input("n_divs_panel", "number of element divisions on panel", 10);
            
            const Real
            length              = _input("panel_l",                                     "length of panel",  0.3),
            tc                  = _input("t_by_c",                                      "length of panel",  0.05),
            ff_to_panel_l       = _input("farfield_to_l_ratio", "Ratio of distance of farfield boundary to panel length",  5.0),
            ff_to_panel_e_size  = _input("farfield_to_panel_elem_size_ratio", "Ratio of element size at far-field to element size at panel",  20.0);

            std::string
            s                   = _input("elem_type",  "type of geometric element in the fluid mesh",     "quad4");
            libMesh::ElemType
            elem_type           = libMesh::Utility::string_to_enum<libMesh::ElemType>(s);

            std::vector<Real>
            x_div_loc           = {-length*ff_to_panel_l, 0., length, length*(1.+ff_to_panel_l)},
            x_relative_dx       = {ff_to_panel_e_size, 1., 1., ff_to_panel_e_size},
            y_div_loc           = {0., length*ff_to_panel_l},
            y_relative_dx       = {1., ff_to_panel_e_size};
            
            std::vector<unsigned int>
            x_divs              = {n_divs_ff_to_panel, n_divs_panel, n_divs_ff_to_panel},
            y_divs              = {n_divs_ff_to_panel};
            
            MAST::MeshInitializer::CoordinateDivisions
            x_coord_divs,
            y_coord_divs;
            
            x_coord_divs.init(nx_divs, x_div_loc, x_relative_dx, x_divs);
            y_coord_divs.init(ny_divs, y_div_loc, y_relative_dx, y_divs);
            
            std::vector<MAST::MeshInitializer::CoordinateDivisions*>
            divs = {&x_coord_divs, &y_coord_divs};
            
            // initialize the mesh
            MAST::PanelMesh2D().init(tc,               // t/c
                                     false,            // if cos bump
                                     1,                // n max bumps
                                     divs,
                                     *_mesh,
                                     elem_type);
        }
        
        if (bc) {
            
            // create the boundary conditions for slip-wall, symmetry and far-field
            MAST::BoundaryConditionBase
            *far_field   = new MAST::BoundaryConditionBase(MAST::FAR_FIELD),
            *symm_wall   = new MAST::BoundaryConditionBase(MAST::SYMMETRY_WALL),
            *slip_wall   = new MAST::BoundaryConditionBase(MAST::SLIP_WALL);
            
            // For Euler flow the panel is modeled with slip wall ...
            _discipline->add_side_load(   4, *slip_wall);
            // ... and the remaining boundary on the bottom of the flow domain is
            // modeled using symmetry wall
            _discipline->add_side_load(   5, *symm_wall);
            // all other boundary (right, top, left) are modeled as far-field
            // conditions.
            for (unsigned int i=1; i<=3; i++)
                _discipline->add_side_load(i, *far_field);

            // store the pointers for later deletion in the destructor
            _boundary_conditions.insert(far_field);
            _boundary_conditions.insert(symm_wall);
            _boundary_conditions.insert(slip_wall);
        }
    }

    void _init_panel_3D(bool mesh, bool bc) {
    
        libmesh_error(); // to be implemented
    }

    //
    // \section fluid_sol_init Initialization of solution
    //
    // \subsection initial_sol Initial solution
    //
    void _init_solution() {
        
        bool
        restart = _input("restart_simulation", "restart simulation from solution", false);
        
        if (!restart) {
            
            std::string
            type = _input("initial_solution", "initial solution field (uniform, random, shock)", "uniform");

            RealVectorX s = RealVectorX::Zero(_dim+2);
            s(0) = _flight_cond->rho();
            s(1) = _flight_cond->rho_u1();
            s(2) = _flight_cond->rho_u2();
            if (_dim > 2)
                s(3) = _flight_cond->rho_u3();
            s(_dim+1) = _flight_cond->rho_e();

            if (type == "uniform")
                // uniform solution
                _sys_init->initialize_solution(s);
            else if (type == "shock") {

                Real
                beta   = _input("beta", "initial shock angle (deg)", 34.799),
                frac   = _input("shock_location_fraction", "location of shock location on panel as a fraction of panel length", 0.5),
                length = _input("panel_l",                                     "length of panel",  0.3);
;
                ShockImpingementSolutions
                *f_sol = new ShockImpingementSolutions(_flight_cond->gas_property.cp,
                                                       _flight_cond->gas_property.cv,
                                                       _flight_cond->mach,
                                                       _flight_cond->gas_property.rho,
                                                       _flight_cond->gas_property.T,
                                                       beta*acos(-1)/180.,
                                                       frac*length);
                _flight_cond->inf_sol.reset(f_sol);
                
                _sys_init->initialize_solution(*f_sol);
            }
            else if (type == "random") {
                
                Real
                amplitude = _input("perturb_amplitude", "maximum amplitude of velocity perturbation as fraction of velocity magnitude", 0.1),
                delta     = amplitude * _flight_cond->velocity_magnitude();
                
                class RandomVelocityPerturbationSolution:
                public MAST::FieldFunction<RealVectorX> {
                    
                protected:
                    
                    unsigned int                 _dim;
                    Real                         _delta;
                    const MAST::FlightCondition &_flt;
                    
                public:
                    
                    RandomVelocityPerturbationSolution(unsigned int dim, Real delta, const MAST::FlightCondition& flt):
                    MAST::FieldFunction<RealVectorX>("sol"), _dim(dim), _delta(delta), _flt(flt) { }
                    
                    virtual ~RandomVelocityPerturbationSolution() { }
                    
                    virtual void
                    operator()(const libMesh::Point& p, const Real t, RealVectorX& v) const {
                        RealVectorX
                        r     = _delta * _flt.rho() * RealVectorX::Random(_dim);
                        v     = RealVectorX::Zero(_dim+2);
                        v(0) = _flt.rho();
                        v(1) = _flt.rho_u1() + r(0);
                        v(2) = _flt.rho_u2() + r(1);
                        if (_dim > 2)
                            v(3) = _flt.rho_u3() + r(2);
                        
                        // updated KE
                        Real
                        ke = 0.;
                        for (unsigned int i=0; i<_dim; i++) ke += std::pow(v(i+1)/v(0), 2);
                        ke = 0.5 * v(0) * std::pow(ke, 0.5);
                        
                        // update the kinetic energy with the random perturbation
                        v(_dim+1) = (_flt.rho_e() - _flt.q0() + ke);
                    }
                };
                
                _sys_init->initialize_solution(RandomVelocityPerturbationSolution(_dim,
                                                                                  delta,
                                                                                  *_flight_cond));
            }
            else
                libmesh_error_msg("Unknown initial solution type: "+type);
        }
        else {
            
            std::string
            output_name = _input("output_file_root", "prefix of output file names", "output"),
            dir_name    = _input("output_file_dir", "directory in which solution vector is stored", "data");
            unsigned int
            t_step      = _input("restart_time_step", "time-step to restart solution", 0);
            
            
            std::ostringstream oss;
            oss << output_name << "_sol_t_" << t_step;
            _sys->read_in_vector(*_sys->solution, dir_name, oss.str(), true);
        }

    }

    
    // \subsection initial_sens_sol Initial sensitivity solution
    void _init_sensitivity_solution() {
        
        bool
        restart = _input("restart_simulation", "restart simulation from solution", false);
        
        if (!restart) {
            
            unsigned int
            n = _mesh->mesh_dimension();
            RealVectorX s = RealVectorX::Zero(n+2);
            /*s(0) = 0.;
            s(1) = _flight_cond->rho_u1_sens_mach();
            s(2) = _flight_cond->rho_u2_sens_mach();
            //if (n > 2)
            //    s(3) = _flight_cond->rho_u3_sens_mach();
            s(n+1) = _flight_cond->rho_e_sens_mach();*/
            s(0) = _flight_cond->rho_sens_rho();
            s(1) = _flight_cond->rho_u1_sens_rho();
            s(2) = _flight_cond->rho_u2_sens_rho();
            //if (n > 2)
            //    s(3) = _flight_cond->rho_u3_sens_mach();
            s(n+1) = _flight_cond->rho_e_sens_rho();
            
            _sys_init->initialize_solution(s);
        }
        else {
            
            _sys->solution->zero();
            _sys->solution->close();
        }
    }

    
public:

    // \section flow_analysis_class Example class
    // \subsection flow_analysis_class_constructor Constructor
    //
    FlowAnalysis(libMesh::LibMeshInit& init,
                 MAST::Examples::GetPotWrapper& input):
    _init           (init),
    _input          (input),
    _dim            (0),
    _mesh           (nullptr),
    _eq_sys         (nullptr),
    _sys            (nullptr),
    _sys_init       (nullptr),
    _discipline     (nullptr),
    _flight_cond    (nullptr),
    _output         (nullptr) {


        _mesh              = new libMesh::ParallelMesh(_init.comm());
        
        // create equation system
        _eq_sys = new libMesh::EquationSystems(*_mesh);
        
        // add the system to be used for fluid analysis
        _sys = &(_eq_sys->add_system<MAST::NonlinearSystem>("fluid"));
        _sys->set_init_B_matrix();
        
        // initialize the mesh. Details of parameters for each mesh are
        // described above.
        _init_mesh(true, false);
        
        // create the discipline where boundary conditions will be stored
        _discipline = new MAST::ConservativeFluidDiscipline(*_eq_sys);
        
        // create system initialization object to add variables
        std::string s;
        s                   = input("fe_order", "order of finite element shape basis functions",     "first");
        libMesh::Order
        fe_order            = libMesh::Utility::string_to_enum<libMesh::Order>(s);
        s                   = input("fe_family",      "family of finite element shape functions", "lagrange");
        libMesh::FEFamily
        fe_family           = libMesh::Utility::string_to_enum<libMesh::FEFamily>(s);

        _sys_init = new MAST::ConservativeFluidSystemInitialization(*_sys,
                                                                    _sys->name(),
                                                                    libMesh::FEType(fe_order, fe_family),
                                                                    _dim);
        
        // set fluid properties
        _flight_cond = new MAST::FlightCondition;
        _flight_cond->flow_unit_vector(0)  =
        _input("flow_unit_vector", "unit vector defining direction of flow", 1., 0);
        _flight_cond->flow_unit_vector(1)  =
        _input("flow_unit_vector", "unit vector defining direction of flow", 0., 1);
        _flight_cond->flow_unit_vector(2)  =
        _input("flow_unit_vector", "unit vector defining direction of flow", 0., 2);
        _flight_cond->mach             = input("mach", "fluid Mach number",                           0.5);
        _flight_cond->gas_property.cp  = input("cp",   "fluid specific heat at constant pressure",  1003.);
        _flight_cond->gas_property.cv  = input("cv",   "fluid specific heat at constant volume",     716.);
        _flight_cond->gas_property.T   = input("T",    "fluid absolute temperature",                 300.);
        _flight_cond->gas_property.rho = input("rho",  "fluid density",                              1.05);
        _flight_cond->gas_property.if_viscous =
        _input("if_viscous", "if the flow analysis should include viscosity", false);
        _flight_cond->init();
        
        // tell the discipline about the fluid values
        _discipline->set_flight_condition(*_flight_cond);

        // initialize the boundary conditions before initialization of the
        // equation system
        _init_mesh(false, true);
        
        // initialize the equation system
        _eq_sys->init();
        
        // print the information
        _mesh->print_info();
        _eq_sys->print_info();

        // initialize the fluid solution
        _init_solution();
        _sys->update();
    }
    
    // \subsection flow_analysis_class_destructor Destructor
    ~FlowAnalysis() {
        
        delete _eq_sys;
        delete _mesh;
        
        delete _discipline;
        delete _sys_init;
        delete _flight_cond;
        
        delete _output;

        {
            std::set<MAST::BoundaryConditionBase*>::iterator
            it   = _boundary_conditions.begin(),
            end  = _boundary_conditions.end();
            for ( ; it!=end; it++)
                delete *it;
        }

        {
            std::set<libMesh::PeriodicBoundary*>::iterator
            it   = _periodic_boundaries.begin(),
            end  = _periodic_boundaries.end();
            for ( ; it!=end; it++)
                delete *it;
        }
    }
    
    // \section flow_computation Computation
    // \subsection flow_transient_analysis  Transient analysis
    void compute_flow() {
        
        bool
        output     = _input("if_output", "if write output to a file", true);
        std::string
        output_name = _input("output_file_root", "prefix of output file names", "output"),
        transient_output_name = output_name + "_transient.exo";
        
        
        // create the nonlinear assembly object
        MAST::TransientAssembly                                  assembly;
        MAST::ConservativeFluidTransientAssemblyElemOperations   elem_ops;
        MAST::FirstOrderNewmarkTransientSolver                   solver;
        RealVectorX
        nvec = RealVectorX::Zero(3);
        nvec(0) =
        _input("force_unit_vector", "unit vector defining direction of integrated force output", 1., 0);
        nvec(1) =
        _input("force_unit_vector", "unit vector defining direction of integrated force output", 0., 1);
        nvec(2) =
        _input("force_unit_vector", "unit vector defining direction of integrated force output", 0., 2);


        MAST::IntegratedForceOutput                              force(nvec);
        std::set<libMesh::boundary_id_type> bids;
        bids.insert(3);
        force.set_participating_boundaries(bids);
        
        assembly.set_discipline_and_system(*_discipline, *_sys_init);
        elem_ops.set_discipline_and_system(*_discipline, *_sys_init);
        force.set_discipline_and_system(*_discipline, *_sys_init);
        solver.set_discipline_and_system(*_discipline, *_sys_init);
        solver.set_elem_operation_object(elem_ops);
        
        //this->initialize_solution();
        
        // file to write the solution for visualization
        libMesh::ExodusII_IO transient_output(*_mesh);
        std::ofstream force_output;
        force_output.open("force.txt");
        force_output
        << std::setw(10) << "t"
        << std::setw(30) << "force" << std::endl;
        
        // time solver parameters
        Real
        factor   = 0.,
        min_fac  = 1.5,
        vel_0    = 0.,
        vel_1    = 1.e12,
        p        = 0.5,
        tval     = 0.,
        max_dt   = _input("max_dt", "maximum time-step size", 1.e-1);
        
        unsigned int
        t_step            = 0,
        iter_count_dt     = 0,
        n_iters_change_dt = _input("n_iters_change_dt", "number of time-steps before dt is changed", 4),
        n_steps           = _input("n_transient_steps", "number of transient time-steps", 100);
        solver.dt         = _input("dt", "time-step size",    1.e-3);
        // the default parameter adds some numerical damping to improve
        // convergence towards steady state solution.
        solver.beta       = _input("beta", "Newmark solver beta parameter ",  0.7);

        // print the fluid values:
        libMesh::out
        << std::setw(15) << "rho: " << std::setw(20) << _flight_cond->rho() << std::endl
        << std::setw(15) << "T: " << std::setw(20) << _flight_cond->gas_property.T << std::endl
        << std::setw(15) << "Mach: " << std::setw(20) << _flight_cond->mach << std::endl
        << std::setw(15) << "gamma: " << std::setw(20) << _flight_cond->gas_property.gamma << std::endl
        << std::setw(15) << "R: " << std::setw(20) << _flight_cond->gas_property.R << std::endl
        << std::setw(15) << "a: " << std::setw(20) << _flight_cond->gas_property.a << std::endl
        << std::setw(15) << "q_dyn: " << std::setw(20) << _flight_cond->q0() << std::endl
        << std::endl;

        // ask the solver to update the initial condition for d2(X)/dt2
        // This is recommended only for the initial time step, since the time
        // integration scheme updates the velocity and acceleration at
        // each subsequent iterate
        solver.solve_highest_derivative_and_advance_time_step(assembly);

        // loop over time steps
        while (t_step < n_steps) {
            
            if (iter_count_dt == n_iters_change_dt) {
                
                libMesh::out
                << "Changing dt:  old dt = " << solver.dt
                << "    new dt = " ;
                
                factor        = std::pow(vel_0/vel_1, p);
                factor        = std::max(factor, min_fac);
                solver.dt     = std::min(solver.dt*factor, max_dt);
                
                libMesh::out << solver.dt << std::endl;
                
                iter_count_dt = 0;
                vel_0         = vel_1;
            }
            else
                iter_count_dt++;
            
            libMesh::out
            << "Time step: "    << t_step
            << " :  t = "       << tval
            << " :  dt = "      << solver.dt
            << " :  xdot-L2 = " << solver.velocity().l2_norm()
            << std::endl;
            
            // write the time-step
            if (output) {
                
                transient_output.write_timestep(transient_output_name,
                                                *_eq_sys,
                                                t_step+1,
                                                _sys->time);
                std::ostringstream oss;
                oss << output_name << "_sol_t_" << t_step;
                _sys->write_out_vector(*_sys->solution, "data", oss.str(), true);
            }
            
            // calculate the output quantity
            force.zero_for_analysis();
            assembly.calculate_output(solver.solution(), true, force);
            force_output
            << std::setw(10) << tval
            << std::setw(30) << force.output_total() << std::endl;
            
            //_sys->adjoint_solve(solver, force, assembly, true);
            //_sys->solution->swap(_sys->get_adjoint_solution(0));
            //adjoint_output.write_timestep("adjoint.exo", *_eq_sys, t_step+1, _sys->time);
            //_sys->solution->swap(_sys->get_adjoint_solution(0));
            
            // solve for the time-step
            solver.solve(assembly);
            solver.advance_time_step();
            
            // update time value
            tval  += solver.dt;
            t_step++;
        }
    }
    
    
    // \subsection flow_transient_sensitivity_analysis  Transient sensitivity analysis
    void
    compute_transient_sensitivity(MAST::Parameter& p) {
        
        bool
        output                = _input("if_output", "if write output to a file", false);
        std::string
        output_name           = _input("output_file_root", "prefix of output file names", "output"),
        transient_output_name = output_name + "_transient_sensitivity_" + p.name() + ".exo",
        nonlinear_sol_root    = output_name+std::string("_sol_t_"),
        nonlinear_sol_dir     = _input("nonlinear_sol_dir", "directory containing the location of nonlinear solutions", "");
        
        // the output from analysis should have been saved for sensitivity
        libmesh_assert(output);
        
        // create the nonlinear assembly object
        MAST::TransientAssembly                                  assembly;
        MAST::ConservativeFluidTransientAssemblyElemOperations   elem_ops;
        MAST::FirstOrderNewmarkTransientSolver                   solver;
        RealVectorX
        nvec = RealVectorX::Zero(3);
        nvec(0) =
        _input("force_unit_vector", "unit vector defining direction of integrated force output", 1., 0);
        nvec(1) =
        _input("force_unit_vector", "unit vector defining direction of integrated force output", 0., 1);
        nvec(2) =
        _input("force_unit_vector", "unit vector defining direction of integrated force output", 0., 2);

        MAST::IntegratedForceOutput                              force(nvec);
        std::set<libMesh::boundary_id_type> bids;
        bids.insert(3);
        force.set_participating_boundaries(bids);
        
        assembly.set_discipline_and_system(*_discipline, *_sys_init);
        elem_ops.set_discipline_and_system(*_discipline, *_sys_init);
        force.set_discipline_and_system(*_discipline, *_sys_init);
        solver.set_discipline_and_system(*_discipline, *_sys_init);
        solver.set_elem_operation_object(elem_ops);
        
        // initialize the solution to zero, or to something that the
        // user may have provided
        _init_solution();
        _sys->update();
        _sys->solution->swap(solver.solution_sensitivity());
        _init_sensitivity_solution();
        _sys->solution->swap(solver.solution_sensitivity());
        
        // file to write the solution for visualization
        libMesh::ExodusII_IO exodus_writer(*_mesh);
        
        // file to write the solution for visualization
        libMesh::ExodusII_IO transient_output(*_mesh);
        std::ofstream force_output;
        force_output.open("force.txt");
        force_output
        << std::setw(10) << "t"
        << std::setw(30) << "force"
        << std::setw(30) << "force_sens" << std::endl;
        
        unsigned int
        t_step            = 0,
        n_steps           = _input("n_transient_steps", "number of transient time-steps", 100);
        solver.dt         = _input("dt", "time-step size",    1.e-3);
        _sys->time        = _input("t_initial", "initial time-step",    0.);
        solver.beta       = _input("beta", "Newmark solver beta parameter ",  0.7);

        // ask the solver to update the initial condition for d2(X)/dt2
        // This is recommended only for the initial time step, since the time
        // integration scheme updates the velocity and acceleration at
        // each subsequent iterate
        solver.solve_highest_derivative_and_advance_time_step(assembly, false);
        solver.solve_highest_derivative_and_advance_time_step_with_sensitivity(assembly, p);
        
        // loop over time steps
        while (t_step < n_steps-1) {
            
            libMesh::out
            << "Time step: " << t_step
            << " :  t = " << _sys->time
            << " :  xdot-L2 = " << solver.velocity().l2_norm()
            << std::endl;
            
            // write the time-step
            if (output) {
                
                _sys->solution->swap(solver.solution_sensitivity());
                exodus_writer.write_timestep(transient_output_name,
                                             *_eq_sys,
                                             t_step+1,
                                             _sys->time);
                _sys->solution->swap(solver.solution_sensitivity());
            }
            
            std::ostringstream oss, oss_sens;
            oss << nonlinear_sol_root << t_step+1;
            _sys->read_in_vector(*_sys->solution, nonlinear_sol_dir, oss.str(), true);
            solver.update_velocity(solver.velocity(), *_sys->solution);
            _sys->update();
            oss_sens << nonlinear_sol_root << t_step << "_sens_t";
            _sys->write_out_vector(/*solver.dt, _sys->time,*/ solver.solution_sensitivity(), "data", oss.str(), true);
            
            // solve for the sensitivity time-step
            force.zero_for_analysis();
            assembly.calculate_output(solver.solution(), true, force);
            assembly.calculate_output_direct_sensitivity(solver.solution(), true,
                                                         &solver.solution_sensitivity(), true,
                                                         p, force);
            force_output
            << std::setw(10) << _sys->time
            << std::setw(30) << force.output_total()
            << std::setw(30) << force.output_sensitivity_total(p) << std::endl;
            
            solver.sensitivity_solve(assembly, p);
            solver.advance_time_step(false);
            solver.advance_time_step_with_sensitivity();
            
            // update time step counter
            t_step++;
        }
        
    }
    

    // \subsection flow_transient_stabilized_sensitivity_analysis  Transient stabilized sensitivity analysis
    void
    compute_transient_stabilized_sensitivity(MAST::Parameter& p) {
        
        bool
        output                = _input("if_output", "if write output to a file", false);
        std::string
        output_name           = _input("output_file_root", "prefix of output file names", "output"),
        transient_output_name = output_name + "_transient_sensitivity_" + p.name() + ".exo",
        nonlinear_sol_root    = output_name+std::string("_sol_t_"),
        nonlinear_sol_dir     = _input("nonlinear_sol_dir", "directory containing the location of nonlinear solutions", "");
        
        // the output from analysis should have been saved for sensitivity
        libmesh_assert(output);
        
        // create the nonlinear assembly object
        MAST::TransientAssembly                                     assembly;
        MAST::ConservativeFluidTransientAssemblyElemOperations      elem_ops;
        MAST::StabilizedFirstOrderNewmarkTransientSensitivitySolver solver;
        RealVectorX
        nvec = RealVectorX::Zero(3);
        nvec(0) =
        _input("force_unit_vector", "unit vector defining direction of integrated force output", 1., 0);
        nvec(1) =
        _input("force_unit_vector", "unit vector defining direction of integrated force output", 0., 1);
        nvec(2) =
        _input("force_unit_vector", "unit vector defining direction of integrated force output", 0., 2);


        MAST::IntegratedForceOutput                              force(nvec);
        std::set<libMesh::boundary_id_type> bids;
        bids.insert(3);
        force.set_participating_boundaries(bids);
        
        assembly.set_discipline_and_system(*_discipline, *_sys_init);
        elem_ops.set_discipline_and_system(*_discipline, *_sys_init);
        force.set_discipline_and_system(*_discipline, *_sys_init);
        solver.set_discipline_and_system(*_discipline, *_sys_init);
        solver.set_elem_operation_object(elem_ops);
        
        _init_solution();
        _sys->update();
        _sys->solution->swap(solver.solution_sensitivity());
        _init_sensitivity_solution();
        _sys->solution->swap(solver.solution_sensitivity());

        // file to write the solution for visualization
        libMesh::ExodusII_IO exodus_writer(*_mesh);
        
        // file to write the solution for visualization
        libMesh::ExodusII_IO transient_output(*_mesh);
        std::ofstream force_output;
        force_output.open("force.txt");
        force_output
        << std::setw(10) << "t"
        << std::setw(30) << "force"
        << std::setw(30) << "force_sens" << std::endl;
        
        unsigned int
        t_step            = _input("initial_transient_step", "first transient step from the transient solution to be used for sensitivity", 0),
        n_steps           = _input("n_transient_steps", "number of transient time-steps", 100);
        solver.dt         = _input("dt", "time-step size",    1.e-3);
        _sys->time        = _input("t_initial", "initial time-step",    0.);
        solver.max_amp    = _input("max_amp", "maximum amplitude for the stabilized sensitivity solver",   1.);
        solver.beta       = _input("sensitivity_beta", "beta for stabilized sensitivity solver",   1.);
        solver.max_index  = n_steps;
        solver.set_nolinear_solution_location(nonlinear_sol_root, nonlinear_sol_dir);
        
        // ask the solver to update the initial condition for d2(X)/dt2
        // This is recommended only for the initial time step, since the time
        // integration scheme updates the velocity and acceleration at
        // each subsequent iterate
        
        // loop over time steps
        while (t_step < n_steps) {
            
            libMesh::out
            << "Time step: " << t_step
            << " :  t = " << _sys->time
            << " :  xdot-L2 = " << solver.velocity().l2_norm()
            << std::endl;
            
            // write the time-step
            if (output) {
                
                _sys->solution->swap(solver.solution_sensitivity());
                exodus_writer.write_timestep(transient_output_name,
                                             *_eq_sys,
                                             t_step+1,
                                             _sys->time);
                _sys->solution->swap(solver.solution_sensitivity());
            }
            
            std::ostringstream oss;
            oss << output_name << "_sol_t_" << t_step;
            oss << "_sens_t";
            _sys->write_out_vector(/*solver.dt, _sys->time,*/ solver.solution_sensitivity(), "data", oss.str(), true);
            
            // solve for the sensitivity time-step
            force.zero_for_analysis();
            assembly.calculate_output(solver.solution(), true, force);
            assembly.calculate_output_direct_sensitivity(solver.solution(), true,
                                                         &solver.solution_sensitivity(), true,
                                                         p, force);
            force_output
            << std::setw(10) << _sys->time
            << std::setw(30) << force.output_total()
            << std::setw(30) << force.output_sensitivity_total(p) << std::endl;
            
            solver.sensitivity_solve(assembly, p);
            
            // update time step counter
            t_step++;
        }
        
    }
    
    
};


// \section flow_driver Main function
int main(int argc, const char** argv)
{
    
    // Initialize libMesh library.
    libMesh::LibMeshInit init(argc, argv);

    // initialize the wrapper to read input parameters. This will check
    // if the executable parameters included a parameter of type
    // `input=${filename}`. If included, then the input parameters will be read
    // from this filename. Otherwise, the parameters will be read from the
    // executable arguments. The wrapper uses default values for parameters
    // if none are provided.
    MAST::Examples::GetPotWrapper
    input(argc, argv, "input");

    bool
    analysis    = input("if_analysis", "whether or not to perform analysis", true),
    sensitivity = input("if_sensitivity", "whether or not to perform sensitivity analysis", false),
    stabilized  = input("if_stabilized_sensitivity", "flag to use standard or stabilized sensitivity analysis", false);

    FlowAnalysis flow(init, input);
    if (analysis)
        flow.compute_flow();
    
    if (sensitivity) {
    MAST::Parameter p("dummy", 0.);
        if (!stabilized)
            flow.compute_transient_sensitivity(p);
        else
            flow.compute_transient_stabilized_sensitivity(p);
    }
    
    // END_TRANSLATE
    return 0;
}
