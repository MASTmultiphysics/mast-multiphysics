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
#include "examples/base/input_wrapper.h"
#include "base/nonlinear_system.h"
#include "base/transient_assembly.h"
#include "base/boundary_condition_base.h"
#include "base/field_function_base.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "base/physics_discipline_base.h"
#include "property_cards/isotropic_material_property_card.h"
#include "property_cards/solid_1d_section_element_property_card.h"
#include "boundary_condition/dirichlet_boundary_condition.h"
#include "heat_conduction/heat_conduction_system_initialization.h"
#include "heat_conduction/heat_conduction_nonlinear_assembly.h"
#include "heat_conduction/heat_conduction_transient_assembly.h"
#include "solver/first_order_newmark_transient_solver.h"

// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_elem_type.h"    // ElemType
#include "libmesh/fe_type.h"           // FEFamily, Order
#include "libmesh/serial_mesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/dof_map.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/nonlinear_solver.h"
#include "libmesh/mesh_generation.h"

// Forward declerations
void compute_transient_solution(MAST::PhysicsDisciplineBase& discipline,
                                MAST::SystemInitialization& sys_init,
                                MAST::Examples::GetPotWrapper& input);

void compute_transient_sensitivity(MAST::PhysicsDisciplineBase& discipline,
                                   MAST::SystemInitialization& sys_init,
                                   MAST::Examples::GetPotWrapper& input,
                                   MAST::Parameter& p);


int main(int argc, const char** argv)
{
    // BEGIN_TRANSLATE Transient conduction of a bar
    //
    //   \tableofcontents
    //
    // This example solves an axial bar extension problem.
    //
    // \section conduction_ex1_init Initialization
    //
    // Initialize libMesh library.
    libMesh::LibMeshInit init(argc, argv);

    MAST::Examples::GetPotWrapper
    input(argc, argv, "input");
    //
    // \subsection conduction_ex1_init_mesh Initialize Mesh and System
    //
    // Create Mesh object on default MPI communicator and generate a line mesh (5 elements, 10 units long).
    //   Note that in libMesh, all meshes are parallel by default in the sense that the equations on the mesh are solved in parallel by PETSc.
    //   A "ReplicatedMesh" is one where all MPI processes have the full mesh in memory, as compared to a "DistributedMesh" where the mesh is
    //   "chunked" up and distributed across processes, each having their own piece.
    libMesh::ReplicatedMesh mesh(init.comm());
    libMesh::MeshTools::Generation::build_line(mesh, 20, 0.0, 10.0);
    mesh.print_info();
    mesh.boundary_info->print_info();
    
    // Create EquationSystems object, which is a container for multiple systems of equations that are defined on a given mesh.
    libMesh::EquationSystems equation_systems(mesh);
    
    // Add system of type MAST::NonlinearSystem (which wraps libMesh::NonlinearImplicitSystem) to the EquationSystems container.
    MAST::NonlinearSystem & system = equation_systems.add_system<MAST::NonlinearSystem>("conduction");
    
    // Create a finite element type for the system. Here we use first order
    // Lagrangian-type finite elements.
    libMesh::FEType fetype(libMesh::FIRST, libMesh::LAGRANGE);
    
    // Initialize the system to the correct set of variables for a conduction
    // analysis. In libMesh this is analogous to adding variables (each with
    // specific finite element type/order to the system for a particular
    // system of equations.
    MAST::HeatConductionSystemInitialization conduction_system(system,
                                                               system.name(),
                                                               fetype);
    
    // Initialize a new conduction discipline using equation_systems.
    MAST::PhysicsDisciplineBase discipline(equation_systems);
        
    // Create and add boundary conditions to the conduction system. This Dirichlet BC fixes
    // the temperature of left end of the bar to a value of 10K. This requires
    // creation of a FieldFunction that provides the value of the constrained
    // variable at given spatial location and time.
    class TempValue: public MAST::FieldFunction<RealVectorX> {
    public:
        TempValue(): MAST::FieldFunction<RealVectorX>("Temp") {}
        virtual ~TempValue() {}
        virtual void operator() (const libMesh::Point& p,
                                 const Real t,
                                 RealVectorX& v) const {
            v.setZero(1);
            v(0) = 10;
        }
    };
    TempValue tval;
    MAST::DirichletBoundaryCondition dirichlet_bc;
    // Conduction system has one variable, Temperature, which is constrained on
    // the left boundary.
    dirichlet_bc.init(0, conduction_system.vars(),
                      // The object takes a pointer to this function to set the
                      // value of temperatue on the left boundary
                      &tval,
                      // tell the object that a system variable order is to be used
                      // and that the sytsem has 1 variable in it.
                      libMesh::SYSTEM_VARIABLE_ORDER, 1);
    // Tell the discipline that the left boundary is constrained by this object.
    discipline.add_dirichlet_bc(0, dirichlet_bc);
    discipline.init_system_dirichlet_bc(system);

    // Initialize the equation system since we now know the size of our
    // system matrices (based on mesh, element type, variables in the
    // conduction_system) as well as the setup of dirichlet boundary conditions.
    // This initialization process is basically a pre-processing step to
    // preallocate storage and spread it across processors.
    equation_systems.init();
    equation_systems.print_info();
    
    //
    // \subsection conduction_ex1_properties Initialize Parameters and Property Cards
    //
    // Create parameters.
    MAST::Parameter zero("zero", 0.0);
    MAST::Parameter    kappa_yy("kappa_yy", 5./6.);
    MAST::Parameter    kappa_zz("kappa_zz", 5./6.);
    MAST::Parameter thickness_y("thy", input("thy", "section y-thickness",  0.06));
    MAST::Parameter thickness_z("thz", input("thz", "section z-thickness",  0.02));
    MAST::Parameter k          ("k",   input("k", "thermal conductivity",   190.));
    MAST::Parameter rho        ("rho", input("rho", "material density",    2700.));
    MAST::Parameter cp         ("cp",  input("cp", "thermal capacitance",   864.));
    MAST::Parameter q          ("flux",input("flux", "heat flux",          -2.e3));

    // Create ConstantFieldFunctions used to spread parameters throughout the model.
    MAST::ConstantFieldFunction   thy_f(        "hy", thickness_y);
    MAST::ConstantFieldFunction   thz_f(        "hz", thickness_z);
    MAST::ConstantFieldFunction hyoff_f(    "hy_off",        zero);
    MAST::ConstantFieldFunction hzoff_f(    "hz_off",        zero);
    MAST::ConstantFieldFunction     k_f(      "k_th",           k);
    MAST::ConstantFieldFunction   rho_f(       "rho",         rho);
    MAST::ConstantFieldFunction    cp_f(        "cp",          cp);
    MAST::ConstantFieldFunction     q_f( "heat_flux",           q);
    MAST::ConstantFieldFunction kappa_yy_f("Kappayy",    kappa_yy);
    MAST::ConstantFieldFunction kappa_zz_f("Kappazz",    kappa_zz);

    // Initialize load.
    MAST::BoundaryConditionBase right_end_flux(MAST::HEAT_FLUX);
    right_end_flux.add(q_f);
    discipline.add_side_load(1, right_end_flux);
    
    // Create the material property card ("card" is NASTRAN lingo) and the
    // relevant parameters to it. An isotropic material needs thermal
    // conductivity (k), density (rho) and thermal capacitance (cp) to
    // describe its behavior.
    MAST::IsotropicMaterialPropertyCard material;
    material.add(k_f);
    material.add(rho_f);
    material.add(cp_f);

    // Create the section property card. Attach all property values.
    MAST::Solid1DSectionElementPropertyCard section;
    section.add(thy_f);
    section.add(thz_f);
    section.add(hyoff_f);
    section.add(hzoff_f);
    section.add(kappa_yy_f);
    section.add(kappa_zz_f);

    // Specify a section orientation point and add it to the section.
    RealVectorX orientation = RealVectorX::Zero(3);
    orientation(1) = 1.0;
    section.y_vector() = orientation;
    
    // Attach material to the card.
    section.set_material(material);
    section.set_diagonal_mass_matrix(true);
    
    // Initialize the section and specify the subdomain in the mesh that it applies to.
    section.init();
    discipline.set_property_for_subdomain(0, section);

    // \section conduction_ex1_computation Computation
    //
    // transient solution
    compute_transient_solution(discipline, conduction_system, input);
    
    // transient solution sensitivity with respect to section y-thickness
    compute_transient_sensitivity(discipline, conduction_system, input, thickness_y);

    return 0;
}


// \subsection conduction_ex1_transient_analysis  Transient analysis
void compute_transient_solution(MAST::PhysicsDisciplineBase& discipline,
                                MAST::SystemInitialization& sys_init,
                                MAST::Examples::GetPotWrapper& input) {
    
    std::string
    output_name = input("output_file_root", "prefix of output file names", "output"),
    transient_output_name = output_name + "_transient.exo";
    
    MAST::NonlinearSystem&       sys = sys_init.system();
    libMesh::EquationSystems& eq_sys = sys.get_equation_systems();
    
    // create the nonlinear assembly object
    MAST::TransientAssembly                                  assembly;
    MAST::HeatConductionTransientAssemblyElemOperations      elem_ops;
    MAST::FirstOrderNewmarkTransientSolver                   solver;
    
    assembly.set_discipline_and_system(discipline, sys_init);
    elem_ops.set_discipline_and_system(discipline, sys_init);
    solver.set_discipline_and_system(discipline, sys_init);
    solver.set_elem_operation_object(elem_ops);

    // initialize the solution to zero, or to something that the
    // user may have provided
    sys.get_dof_map().enforce_constraints_exactly(sys, sys.solution.get());
    sys.update();

    // file to write the solution for visualization
    libMesh::ExodusII_IO transient_output(sys.get_mesh());

    unsigned int
    t_step            = 0,
    n_steps           = input("n_transient_steps", "number of transient time-steps", 100);
    solver.dt         = input("dt", "time-step size",    1.e+3);
    solver.beta       = input("beta", "Newmark solver beta parameter ",  0.5);
        
    // ask the solver to update the initial condition for d2(X)/dt2
    // This is recommended only for the initial time step, since the time
    // integration scheme updates the velocity and acceleration at
    // each subsequent iterate
    solver.solve_highest_derivative_and_advance_time_step(assembly);
    
    // loop over time steps
    while (t_step < n_steps) {
        
        libMesh::out
        << "Time step: "    << t_step
        << " :  t = "       << sys.time
        << " :  dt = "      << solver.dt
        << " :  xdot-L2 = " << solver.velocity().l2_norm()
        << std::endl;
        
        // write the time-step
        transient_output.write_timestep(transient_output_name,
                                        eq_sys,
                                        t_step+1,
                                        sys.time);
        std::ostringstream oss_sol;
        oss_sol << output_name << "_sol_t_" << t_step;
        sys_init.system().write_out_vector(*sys.solution, "data", oss_sol.str(), true);
        
        // solve for the time-step
        solver.solve(assembly);
        solver.advance_time_step();

        // update time step
        t_step++;
    }
}


// \subsection conduction_ex1_transient_sensitivity_analysis  Transient sensitivity analysis
void compute_transient_sensitivity(MAST::PhysicsDisciplineBase& discipline,
                                   MAST::SystemInitialization& sys_init,
                                   MAST::Examples::GetPotWrapper& input,
                                   MAST::Parameter& p) {
    
    std::string
    output_name           = input("output_file_root", "prefix of output file names", "output"),
    transient_output_name = output_name + "_transient_sensitivity_" + p.name() + ".exo",
    nonlinear_sol_dir     = input("nonlinear_sol_dir", "directory containing the location of nonlinear solutions", "data");
    
    MAST::NonlinearSystem&       sys = sys_init.system();
    libMesh::EquationSystems& eq_sys = sys.get_equation_systems();

    // create the nonlinear assembly object
    MAST::TransientAssembly                                  assembly;
    MAST::HeatConductionTransientAssemblyElemOperations      elem_ops;
    MAST::FirstOrderNewmarkTransientSolver                   solver;
    
    assembly.set_discipline_and_system(discipline, sys_init);
    elem_ops.set_discipline_and_system(discipline, sys_init);
    solver.set_discipline_and_system(discipline, sys_init);
    solver.set_elem_operation_object(elem_ops);
    
    // initial condition for solution is read from disk and initial condition
    // for sensitivity of temperature is also assumed to be zero.
    std::ostringstream oss;
    oss << output_name << "_sol_t_" << 0;
    sys.read_in_vector(*sys.solution, nonlinear_sol_dir, oss.str(), true);
    sys.update();
    
    // file to write the solution for visualization
    libMesh::ExodusII_IO transient_output(sys.get_mesh());
        
    unsigned int
    t_step            = 0,
    n_steps           = input("n_transient_steps", "number of transient time-steps", 100);
    solver.dt         = input("dt", "time-step size",    1.e+3);
    sys.time          = 0.;
    solver.beta       = input("beta", "Newmark solver beta parameter ",  0.5);
    
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
        << " :  t = " << sys.time
        << " :  xdot-L2 = " << solver.velocity_sensitivity().l2_norm()
        << std::endl;
        
        // write the time-step
        sys.solution->swap(solver.solution_sensitivity());
        transient_output.write_timestep(transient_output_name,
                                        eq_sys,
                                        t_step+1,
                                        sys.time);
        sys.solution->swap(solver.solution_sensitivity());
        
        std::ostringstream oss_sol;
        oss_sol << output_name << "_sol_t_" << t_step+1;
        sys.read_in_vector(*sys.solution, nonlinear_sol_dir, oss_sol.str(), true);
        solver.update_velocity(solver.velocity(), *sys.solution);
        sys.update();


        // solve for the sensitivity time-step
        solver.sensitivity_solve(assembly, p);
        solver.advance_time_step(false);
        solver.advance_time_step_with_sensitivity();

        // update time step counter
        t_step++;
    }
    
}

// END_TRANSLATE
