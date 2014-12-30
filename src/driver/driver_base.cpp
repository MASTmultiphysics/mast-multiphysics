
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


