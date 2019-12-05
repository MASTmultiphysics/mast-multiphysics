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

// libMesh includes.
#include <libmesh/libmesh.h>
#include <libmesh/parallel.h>
#include <libmesh/replicated_mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/equation_systems.h>
#include <libmesh/fe_type.h>
#include <libmesh/numeric_vector.h>

// MAST includes.
#include "base/nonlinear_system.h"
#include "elasticity/structural_system_initialization.h"
#include "base/physics_discipline_base.h"
#include "boundary_condition/dirichlet_boundary_condition.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "property_cards/isotropic_material_property_card.h"
#include "property_cards/solid_1d_section_element_property_card.h"
#include "base/nonlinear_implicit_assembly.h"
#include "elasticity/structural_nonlinear_assembly.h"

int main(int argc, const char** argv)
{
    // BEGIN_TRANSLATE Extension of bar
    //
    // This example solves an axial bar extension problem. 
    //
    // Initialize libMesh library.
    libMesh::LibMeshInit init(argc, argv);

    // Create Mesh object on default MPI communicator and generate a line mesh (5 elements, 10 units long).
    //   Note that in libMesh, all meshes are parallel by default in the sense that the equations on the mesh are solved in parallel by PETSc.
    //   A "ReplicatedMesh" is one where all MPI processes have the full mesh in memory, as compared to a "DistributedMesh" where the mesh is
    //   "chunked" up and distributed across processes, each having their own piece.
    libMesh::ReplicatedMesh mesh(init.comm());
    libMesh::MeshTools::Generation::build_line(mesh, 5, 0.0, 10.0);
    mesh.print_info();
    mesh.boundary_info->print_info();

    // Create EquationSystems object, which is a container for multiple systems of equations that are defined on a given mesh.
    libMesh::EquationSystems equation_systems(mesh);

    // Add system of type MAST::NonlinearSystem (which wraps libMesh::NonlinearImplicitSystem) to the EquationSystems container.
    //   We name the system "structural" and also get a reference to the system so we can easily reference it later.
    MAST::NonlinearSystem & system = equation_systems.add_system<MAST::NonlinearSystem>("structural");

    // Create a finite element type for the system. Here we use first order
    // Lagrangian-type finite elements.
    libMesh::FEType fetype(libMesh::FIRST, libMesh::LAGRANGE);

    // Initialize the system to the correct set of variables for a structural
    // analysis. In libMesh this is analogous to adding variables (each with
    // specific finite element type/order to the system for a particular
    // system of equations.
    MAST::StructuralSystemInitialization structural_system(system,
                                                           system.name(),
                                                           fetype);

    // Initialize a new structural discipline using equation_systems.
    MAST::PhysicsDisciplineBase discipline(equation_systems);

    // Create and add boundary conditions to the structural system. A Dirichlet BC fixes the left end of the bar.
    // This definition uses the numbering created by the libMesh mesh generation function.
    MAST::DirichletBoundaryCondition dirichlet_bc;
    dirichlet_bc.init(0, structural_system.vars());
    discipline.add_dirichlet_bc(0, dirichlet_bc);
    discipline.init_system_dirichlet_bc(system);

    // Initialize the equation system since we now know the size of our
    // system matrices (based on mesh, element type, variables in the
    // structural_system) as well as the setup of dirichlet boundary conditions.
    // This initialization process is basically a pre-processing step to
    // preallocate storage and spread it across processors.
    equation_systems.init();
    equation_systems.print_info();

    // Create parameters.
    MAST::Parameter thickness_y("thy", 0.06);
    MAST::Parameter thickness_z("thz", 0.02);
    MAST::Parameter E("E", 72.0e9);
    MAST::Parameter nu("nu", 0.33);
    MAST::Parameter zero("zero", 0.0);
    MAST::Parameter pressure("p", 2.0e4);

    // Create ConstantFieldFunctions used to spread parameters throughout the model.
    MAST::ConstantFieldFunction thy_f("hy", thickness_y);
    MAST::ConstantFieldFunction thz_f("hz", thickness_z);
    MAST::ConstantFieldFunction E_f("E", E);
    MAST::ConstantFieldFunction nu_f("nu", nu);
    MAST::ConstantFieldFunction hyoff_f("hy_off", zero);
    MAST::ConstantFieldFunction hzoff_f("hz_off", zero);
    MAST::ConstantFieldFunction pressure_f("pressure", pressure);

    // Initialize load.
    // TODO - Switch this to a concentrated/point load on the right end of the bar.
    MAST::BoundaryConditionBase right_end_pressure(MAST::SURFACE_PRESSURE);
    right_end_pressure.add(pressure_f);
    discipline.add_side_load(1, right_end_pressure);

    // Create the material property card ("card" is NASTRAN lingo) and the relevant parameters to it. An isotropic
    // material needs elastic modulus (E) and Poisson ratio (nu) to describe its behavior.
    MAST::IsotropicMaterialPropertyCard material;
    material.add(E_f);
    material.add(nu_f);

    // Create the section property card. Attach all property values.
    MAST::Solid1DSectionElementPropertyCard section;
    section.add(thy_f);
    section.add(thz_f);
    section.add(hyoff_f);
    section.add(hzoff_f);

    // Specify a section orientation point and add it to the section.
    RealVectorX orientation = RealVectorX::Zero(3);
    orientation(1) = 1.0;
    section.y_vector() = orientation;

    // Attach material to the card.
    section.set_material(material);

    // Initialize the section and specify the subdomain in the mesh that it applies to.
    section.init();
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

    // Solve the system and print displacement degrees-of-freedom to screen.
    nonlinear_system.solve(elem_ops, assembly);
    nonlinear_system.solution->print_global();

    // END_TRANSLATE
    return 0;
}
