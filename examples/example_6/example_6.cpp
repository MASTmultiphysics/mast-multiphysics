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
#include "fluid/conservative_fluid_system_initialization.h"
#include "fluid/conservative_fluid_discipline.h"
#include "fluid/conservative_fluid_transient_assembly.h"
#include "fluid/flight_condition.h"
#include "fluid/integrated_force_output.h"
#include "solver/first_order_newmark_transient_solver.h"

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


//
// BEGIN_TRANSLATE Flow analysis
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

    // This will initialize the mesh
    void _init_mesh(bool mesh, bool bc) {
        
        if (mesh) {
            
            libmesh_assert(!_mesh);
            
            _mesh              = new libMesh::ParallelMesh(_init.comm());
        }

        // The mesh is created using classes written in MAST. The particular
        // mesh to be used can be selected using the input parameter
        // ` mesh=val `, where `val` can be one of the following:
        //   - `naca0012` for flow over a NACA 0012 airfoil
        //   - `cylinder` for flow over a cylinder
        //   - `naca0012_wing` for flow over swept wing with NACA0012 section
        //   - `panel_2D` for 2D flow analysis over a panel
        //   - `panel_3D` for 3D flow analysis over a panel
        //
        // The meshing and boundary conditions for each flow analysis case
        // can be modified using suitable parameters, which are described in this
        // example.
        std::string
        s  = _input("mesh",
                    "type of mesh to be analyzed {naca0012, cylinder, naca0012_wing, panel_2D, panel_3D}",
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
        else
            libmesh_error_msg("unknown mesh type");
    }
    
    // If \p mesh is \p true then the mesh will be initialized. If \p bc
    // is true then the boundary conditions will be initialized
    void _init_naca0012(bool mesh, bool bc) {
        
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
                                                  e_type);

        }
        
        if (bc) {

            libmesh_assert(_flight_cond);
            
            // this assumes a viscous analysis
            libmesh_assert(_flight_cond->gas_property.if_viscous);
            
            // create the boundary conditions for slip-wall, symmetry and far-field
            MAST::BoundaryConditionBase
            *far_field   = new MAST::BoundaryConditionBase(MAST::FAR_FIELD),
            *wall        = new MAST::BoundaryConditionBase(MAST::NO_SLIP_WALL);
            
            // Cylinder surface is boundary id 3 ...
            _discipline->add_side_load(   3, *wall);
            // boundary id 1 is far field
            _discipline->add_side_load(   1, *far_field);

            // store the pointers for later deletion in the destructor
            _boundary_conditions.insert(far_field);
            _boundary_conditions.insert(wall);

        }
    }

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
            MAST::PanelMesh2D().init(0.,               // t/c
                                     false,            // if cos bump
                                     0,                // n max bumps
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
    
    void _init_solution() {
        
        bool
        restart = _input("restart_simulation", "restart simulation from solution", false);
        
        if (!restart) {
            
            unsigned int
            n = _mesh->mesh_dimension();
            RealVectorX s = RealVectorX::Zero(n+2);
            s(0) = _flight_cond->rho();
            s(1) = _flight_cond->rho_u1();
            s(2) = _flight_cond->rho_u2();
            if (n > 2)
                s(3) = _flight_cond->rho_u3();
            s(n+1) = _flight_cond->rho_e();
            
            _sys_init->initialize_solution(s);
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

public:

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
    _output         (nullptr) {


        // initialize the mesh. Details of parameters for each mesh are
        // described above.
        _init_mesh(true, false);
        
        // create equation system
        _eq_sys = new libMesh::EquationSystems(*_mesh);
        
        // add the system to be used for fluid analysis
        _sys = &(_eq_sys->add_system<MAST::NonlinearSystem>("fluid"));
        
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
        _flight_cond->gas_property.rho = input("rho",  "fluid density",                              1.35);
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
    }
    
    
    ~FlowAnalysis() {
        
        delete _eq_sys;
        delete _mesh;
        
        delete _discipline;
        delete _sys_init;
        delete _flight_cond;
        
        delete _output;

        std::set<MAST::BoundaryConditionBase*>::iterator
        it   = _boundary_conditions.begin(),
        end  = _boundary_conditions.end();
        for ( ; it!=end; it++)
            delete *it;
    }
    
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
        nvec(1) = 1.;
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
        libMesh::out << "q_dyn = " << _flight_cond->q0() << std::endl;
        
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
            assembly.calculate_output(solver.solution(), force);
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
};


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

    FlowAnalysis flow(init, input);
    flow.compute_flow();
    
    // END_TRANSLATE
    return 0;
}
