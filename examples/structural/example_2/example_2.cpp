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
#include <fstream>

// libMesh includes.
#include <libmesh/libmesh.h>
#include <libmesh/parallel.h>
#include <libmesh/replicated_mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/equation_systems.h>
#include <libmesh/fe_type.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/point_locator_base.h>

// MAST includes.
#include "examples/base/input_wrapper.h"
#include "base/nonlinear_system.h"
#include "elasticity/structural_system_initialization.h"
#include "base/physics_discipline_base.h"
#include "boundary_condition/dirichlet_boundary_condition.h"
#include "boundary_condition/point_load_condition.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "property_cards/isotropic_material_property_card.h"
#include "property_cards/solid_2d_section_element_property_card.h"
#include "base/nonlinear_implicit_assembly.h"
#include "elasticity/structural_nonlinear_assembly.h"
#include "solver/arclength_continuation_solver.h"
#include "solver/pseudo_arclength_continuation_solver.h"

class PointLoad: public MAST::FieldFunction<RealVectorX> {
public:
    PointLoad(MAST::Parameter& p):
    MAST::FieldFunction<RealVectorX>("load"), _p(p) {}
    virtual void operator()(const libMesh::Point& p, const Real t, RealVectorX& v) const {
        v = RealVectorX::Zero(6);
        // force in z direction
        v(2) = _p();
    }
    virtual void derivative(const MAST::FunctionBase& f,
                            const libMesh::Point& p, const Real t, RealVectorX& v) const {
        v = RealVectorX::Zero(6);
        if (&f == &_p)
            v(2) = 1.;
    }

protected:
    MAST::Parameter& _p;
};

int main(int argc, char* argv[])
{
    // BEGIN_TRANSLATE Nonlinear plate bending with continuation solver
    //
    // This example computes the nonlinaer deformation of a curved plate
    // with a temperature and pressure load using path continuation solver.
    //
    // For direct solver, run with options: -ksp_type preonly -pc_type lu
    //
    // Initialize libMesh library.
    libMesh::LibMeshInit init(argc, argv);

    // wapper to GetPot used to read in values for parameters
    MAST::Examples::GetPotWrapper
    input(argc, argv, "input");

    
    // Create Mesh object on default MPI communicator and generate a 2D mesh of QUAD4 elements. We discretize with
    // 12 elements in the x-direction (0.0 to 36.0 inches) and 40 elements in the y-direction (0.0 to 120.0 inches).
    //   Note that in libMesh, all meshes are parallel by default in the sense that the equations on the mesh are
    //   solved in parallel by PETSc. A "ReplicatedMesh" is one where all MPI processes have the full mesh in memory,
    //   as compared to a "DistributedMesh" where the mesh is "chunked" up and distributed across processes, each
    //   having their own piece.
    Real
    length = input("length", "length of the plate" ,.508),
    width  = input("width", "length of the plate"  ,.508),
    R      = input("radius", "radius of the circle where the circumference defines the curved plate", 2.540);
    
    unsigned int
    nx     = input("nx", "number of elements along the x-axis" , 20),
    ny     = input("ny", "number of elements along the y-axis" , 20);

    libMesh::ReplicatedMesh mesh(init.comm());
    libMesh::MeshTools::Generation::build_square(mesh, nx, ny, 0.0, length, 0.0, width, libMesh::QUAD9);
    mesh.print_info();

    // get the node for which the displacement will be written to the
    // load.txt file.
    libMesh::Point
    pt(length/2., width/2., 0.), // location of mid-point before shift
    pt0,
    dr1, dr2;
    const libMesh::Node
    *nd = nullptr;

    // if a finite radius is defined, change the mesh to a circular arc of specified radius
    libMesh::MeshBase::node_iterator
    n_it   = mesh.nodes_begin(),
    n_end  = mesh.nodes_end();

    // initialize the pointer to a node
    nd   = *n_it;
    pt0  = *nd;
    
    for ( ; n_it != n_end; n_it++) {
        
        dr1  = pt0;
        dr1 -= pt;
        
        dr2  = **n_it;
        dr2 -= pt;
        
        if (dr2.norm() < dr1.norm()) {

            nd  = *n_it;
            pt0 = *nd;
        }
        
        // compute angle based on y-location
        if (R > 0.) {
            Real a  = asin(((**n_it)(1)-width*.5)/R);
            (**n_it)(2) = cos(a)*R;
        }
    }
    
    std::cout << *nd << std::endl;
    
    // also, create a boundary id for this node since we will constrain the
    // x-displacement here
    mesh.boundary_info->add_node(nd, 6);
    mesh.boundary_info->sideset_name(6) = "mid_node";
    mesh.prepare_for_use();

    // For later reference, the boundary and subdomain ID's generated by the libMesh mesh generation are sketched below.
    //
    // \verbatim
    //        y               (#) Boundary ID
    //        | (2)           [#] Subdomain ID
    //  120.0 O-----O
    //        |     |          x - right
    //        |     |          y - up
    //    (3) | [0] | (1)      z - out of screen
    //        |     |
    //        |     |
    //    0.0 O-----O--x
    //          (0)
    //       0.0   36.0
    // \endverbatim
    //
    // Create EquationSystems object, which is a container for multiple systems of equations that are
    // defined on a given mesh.
    libMesh::EquationSystems equation_systems(mesh);

    // Add system of type MAST::NonlinearSystem (which wraps libMesh::NonlinearImplicitSystem) to the
    // EquationSystems container.
    //   We name the system "structural" and also get a reference to the system so we can easily reference it later.
    MAST::NonlinearSystem & system = equation_systems.add_system<MAST::NonlinearSystem>("structural");

    // Create a finite element type for the system. Here we use first order Lagrangian-type finite elements.
    libMesh::FEType fetype(libMesh::SECOND, libMesh::LAGRANGE);

    // Initialize the system to the correct set of variables for a structural analysis. In libMesh this is analogous
    // to adding variables (each with specific finite element type/order to the system for a particular system
    // of equations.
    MAST::StructuralSystemInitialization structural_system(system,
                                                           system.name(),
                                                           fetype);

    // Initialize a new structural discipline using equation_systems.
    MAST::PhysicsDisciplineBase discipline(equation_systems);

    // Create and add boundary conditions to the structural system. A Dirichlet
    // BC adds fixed displacement BCs. Here we
    // use the side boundary ID numbering created by the libMesh generator to
    // clamp the edge of the mesh along x=0.0.
    // We apply the simply supported BC on bottom and top. For clamped edges,
    // set vars={0, 1, 2, 3, 4, 5, 6}.
    MAST::DirichletBoundaryCondition clamped_edge0, clamped_edge1, clamped_edge2, clamped_edge3;    // Create BC object.
    std::vector<unsigned int>
    vars   = {1, 2},
    vars_x = {0};
    clamped_edge0.init(0, vars);   // Assign boundary ID and variables to constrain
    clamped_edge1.init(6, vars_x);   // Assign boundary ID and variables to constrain
    clamped_edge2.init(2, vars);   // Assign boundary ID and variables to constrain
    clamped_edge3.init(3, vars);   // Assign boundary ID and variables to constrain
    
    // this applies the simply supported condition to the bottom (bid=0) and
    // top (bid=1) edges. For setting conditions to left (bid=3) and right (bid=1)
    // call the method with the respective boundary ids.
    discipline.add_dirichlet_bc(0, clamped_edge0);     // Attach boundary condition to discipline
    discipline.add_dirichlet_bc(6, clamped_edge1);     // Attach boundary condition to discipline
    discipline.add_dirichlet_bc(2, clamped_edge2);     // Attach boundary condition to discipline
    //discipline.add_dirichlet_bc(3, clamped_edge3);     // Attach boundary condition to discipline
    discipline.init_system_dirichlet_bc(system);      // Initialize the BC in the system.

    // Initialize the equation system since we now know the size of our system matrices (based on mesh, element type,
    // variables in the structural_system) as well as the setup of dirichlet boundary conditions. This initialization
    // process is basically a pre-processing step to preallocate storage and spread it across processors.
    equation_systems.init();
    equation_systems.print_info();

    // Create parameters.
    MAST::Parameter thickness  ("th",    input(   "th",      "plate thickness",              .0127));
    MAST::Parameter E          ("E",     input(    "E",      "Young's modulus",          3.10275e9));
    MAST::Parameter nu         ("nu",    input(   "nu",     "Poission's ratio",                 .3));
    MAST::Parameter rho        ("rho",   input(  "rho",     "material density",        0.1*0.00259));
    MAST::Parameter kappa      ("kappa", input("kappa",          "shear correction factor",  5./6.));
    MAST::Parameter alpha      ("alpha", input("alpha", "coefficient of thermal expansion", 1.5e-5));
    MAST::Parameter zero       ("zero", 0.0);
    MAST::Parameter pressure   ("p",     input(    "p",     "initial point load",                0));
    MAST::Parameter temperature("T",     input( "temp",  "initial temperature",                0.0));

    // Create ConstantFieldFunctions used to spread parameters throughout the model.
    MAST::ConstantFieldFunction th_f("h", thickness);
    MAST::ConstantFieldFunction E_f("E", E);
    MAST::ConstantFieldFunction nu_f("nu", nu);
    MAST::ConstantFieldFunction rho_f("rho", rho);
    MAST::ConstantFieldFunction kappa_f("kappa", kappa);
    MAST::ConstantFieldFunction alpha_f("alpha_expansion", alpha);
    MAST::ConstantFieldFunction off_f("off", zero);
    MAST::ConstantFieldFunction temperature_f("temperature", temperature);
    MAST::ConstantFieldFunction ref_temp_f("ref_temperature", zero);
    PointLoad load_f(pressure);
    
    // Initialize point load.
    MAST::PointLoadCondition point_load(MAST::POINT_LOAD);
    point_load.add(load_f);
    point_load.add_node(*nd);
    discipline.add_point_load(point_load);

    MAST::BoundaryConditionBase temperature_load(MAST::TEMPERATURE);
    temperature_load.add(temperature_f);
    temperature_load.add(ref_temp_f);
    discipline.add_volume_load(0, temperature_load);

    // Create the material property card ("card" is NASTRAN lingo) and the relevant parameters to it. An isotropic
    // material needs elastic modulus (E) and Poisson ratio (nu) to describe its behavior.
    MAST::IsotropicMaterialPropertyCard material;
    material.add(E_f);
    material.add(nu_f);
    material.add(alpha_f);
    material.add(rho_f);

    // Create the section property card. Attach all property values.
    MAST::Solid2DSectionElementPropertyCard section;
    section.add(th_f);
    section.add(off_f);
    section.add(kappa_f);
    section.set_strain(MAST::NONLINEAR_STRAIN);

    // Attach material to the card.
    section.set_material(material);

    // Initialize the specify the subdomain in the mesh that it applies to.
    discipline.set_property_for_subdomain(0, section);

    // Create nonlinear assembly object and set the discipline and
    // structural_system. Create reference to system.
    MAST::NonlinearImplicitAssembly assembly;
    MAST::StructuralNonlinearAssemblyElemOperations elem_ops;
    assembly.set_discipline_and_system(discipline, structural_system);
    elem_ops.set_discipline_and_system(discipline, structural_system);
    MAST::NonlinearSystem& nonlinear_system = assembly.system();

    // Zero the solution before solving.
    nonlinear_system.solution->zero();

    // Write output to Exodus.
    libMesh::ExodusII_IO exodus_io(mesh);
    
    const unsigned int
    dof_num = nd->dof_number(0, 2, 0);
    
    // Solve the system and print displacement degrees-of-freedom to screen.
    unsigned int
    n_temp_steps  = input( "n_temp_steps", "number of load steps for temperature increase",  10),
    n_press_steps = input("n_press_steps",    "number of load steps for pressure increase",  30);

    // write the header to the load.txt file
    Real
    max_temp = input("max_temp", "maximum temperature", 0.),
    dt = 1./(n_temp_steps+n_press_steps-1.);
    std::ofstream out;
    if (mesh.comm().rank() == 0) {
        out.open("load.txt", std::ofstream::out);
        out
        << std::setw(10) << "iter"
        << std::setw(25) << "temperature"
        << std::setw(25) << "pressure"
        << std::setw(25) << "displ" << std::endl;
    }
    
    // first solve the the temperature increments
    std::vector<Real> vec1;
    std::vector<unsigned int> vec2 = {dof_num};
    if (n_temp_steps) {
        
        MAST::PseudoArclengthContinuationSolver solver;
        solver.schur_factorization = input("if_schur_factorization", "use Schur-factorization in continuation solver", true);
        solver.min_step            = input("min_step", "minimum arc-length step-size for continuation solver",          10.);

        // specify temperature as the load parameter to be changed per
        // load step
        solver.set_assembly_and_load_parameter(elem_ops, assembly, temperature);
        
        // the initial deformation direction is identified with a
        // unit change in temperature.
        solver.initialize(temperature());
        // with the search direction defined, we define the arc length
        // per load step to be a factor of 2 greater than the initial step.
        solver.arc_length *= 2;
        
        for (unsigned int i=0; i<n_temp_steps; i++) {

            solver.solve();
            libMesh::out
            << "  iter: " << i
            << "  temperature: " << temperature()
            << "  pressure: "    << pressure() << std::endl;
            
            // get the value of the node at the center of the plate for output
            system.solution->localize(vec1, vec2);
            
            // write the value to the load.txt file
            if (mesh.comm().rank() == 0) {
                out
                << std::setw(10) << i
                << std::setw(25) << temperature()
                << std::setw(25) << pressure()
                << std::setw(25) << vec1[0] << std::endl;
            }
            system.time += dt;
            
            // write the current solution to the exodus file for
            // visualization
            exodus_io.write_timestep("ex3_plate_static.exo",
                                     equation_systems,
                                     i+1,
                                     system.time);
            
            if (temperature() > max_temp) {
                temperature = max_temp;
                nonlinear_system.solve(elem_ops, assembly);
                break;
            }
        }
    }

    
    // next, solve the the pressure increments
    if (n_press_steps) {
        
        MAST::PseudoArclengthContinuationSolver solver;
        solver.schur_factorization = input("if_schur_factorization", "use Schur-factorization in continuation solver", true);
        solver.min_step            = input("min_step", "minimum arc-length step-size for continuation solver",          10.);

        // specify pressure as the load parameter to be changed per
        // load step
        solver.set_assembly_and_load_parameter(elem_ops, assembly, pressure);
        
        // initial search direction is defined usign a pressure of 2e3.
        solver.initialize(-1.e1);
        
        // the arch length is chanegd to a factor of 4 for each load step.
        solver.arc_length *= 10;

        for (unsigned int i=0+n_temp_steps; i<n_press_steps+n_temp_steps; i++) {
            solver.solve();
            std::cout
            << "  iter: " << i
            << "  temperature: " << temperature()
            << "  pressure: "    << pressure() << std::endl;
            system.solution->localize(vec1, vec2);
            if (mesh.comm().rank() == 0) {
                out
                << std::setw(10) << i
                << std::setw(25) << temperature()
                << std::setw(25) << pressure()
                << std::setw(25) << vec1[0] << std::endl;
            }
            system.time += dt;
            exodus_io.write_timestep("ex3_plate_static.exo",
                                     equation_systems,
                                     i+1,
                                     system.time);
        }
    }

    // END_TRANSLATE
    return 0;
}
