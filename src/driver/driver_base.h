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

#ifndef __mast__driver_base__
#define __mast__driver_base__

namespace MAST {
    // Forward declerations
    class PhysicsDisciplineBase;
    class SystemInitialization;
    class NonlinearImplicitAssembly;
    class TransientAssembly;
    class TransientSolverBase;
    class OutputAssemblyBase;
    class Parameter;
    

    namespace Driver {
        
        bool nonlinear_solution(MAST::PhysicsDisciplineBase&     discipline,
                                MAST::SystemInitialization&      system,
                                MAST::NonlinearImplicitAssembly& assembly);


        bool sensitivity_solution(MAST::PhysicsDisciplineBase&     discipline,
                                  MAST::SystemInitialization&      system,
                                  MAST::NonlinearImplicitAssembly& assembly,
                                  MAST::Parameter&        f);

        
        bool transient_solution_step(MAST::PhysicsDisciplineBase&     discipline,
                                     MAST::SystemInitialization&      system,
                                     MAST::TransientAssembly&         assembly,
                                     MAST::TransientSolverBase&       solver);

        bool output_evaluation(MAST::PhysicsDisciplineBase&     discipline,
                               MAST::SystemInitialization&      system,
                               MAST::OutputAssemblyBase&        assembly);

    
        bool adjoint_solution(MAST::PhysicsDisciplineBase&     discipline,
                              MAST::SystemInitialization&      system,
                              MAST::NonlinearImplicitAssembly& nonlin_assembly,
                              MAST::OutputAssemblyBase&        output_assembly);
    }
}


#endif // __mast__driver_base__
