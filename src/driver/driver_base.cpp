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

// MAST includes
#include "driver/driver_base.h"
#include "base/nonlinear_implicit_assembly.h"
#include "base/transient_assembly.h"
#include "solver/transient_solver_base.h"
#include "base/output_assembly_base.h"
#include "base/physics_discipline_base.h"


// libMesh includes
#include "libmesh/nonlinear_implicit_system.h"


bool
MAST::Driver::nonlinear_solution(MAST::PhysicsDisciplineBase&     discipline,
                                 MAST::SystemInitialization&      system,
                                 MAST::NonlinearImplicitAssembly& assembly) {
    
    assembly.attach_discipline_and_system(discipline, system);
    
    libMesh::NonlinearImplicitSystem&      nonlin_sys   =
    dynamic_cast<libMesh::NonlinearImplicitSystem&>(assembly.system());
    
    nonlin_sys.solve();
    
    assembly.clear_discipline_and_system();
    
    return true;
}



bool
MAST::Driver::transient_solution_step(MAST::PhysicsDisciplineBase&     discipline,
                                      MAST::SystemInitialization&      system,
                                      MAST::TransientAssembly&         assembly,
                                      MAST::TransientSolverBase&       solver) {
    
    assembly.attach_discipline_and_system(discipline, solver, system);
    solver.set_assembly(assembly);
    
    solver.solve();

    assembly.clear_discipline_and_system();
    
    return true;
}




bool
MAST::Driver::output_evaluation(MAST::PhysicsDisciplineBase&     discipline,
                                MAST::SystemInitialization&      system,
                                MAST::OutputAssemblyBase&        assembly) {
    
    assembly.attach_discipline_and_system(discipline, system);

    assembly.system().assemble_qoi();
    
    assembly.clear_discipline_and_system();
    
    return true;
}



bool
MAST::Driver::adjoint_solution(MAST::PhysicsDisciplineBase&     discipline,
                               MAST::SystemInitialization&      system,
                               MAST::NonlinearImplicitAssembly& nonlin_assembly,
                               MAST::OutputAssemblyBase&        output_assembly) {
    
    nonlin_assembly.attach_discipline_and_system(discipline, system);
    output_assembly.attach_discipline_and_system(discipline, system);
    
    output_assembly.system().adjoint_solve();
    
    nonlin_assembly.clear_discipline_and_system();
    output_assembly.clear_discipline_and_system();
    
    return true;
}


