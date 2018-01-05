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

#ifndef __mast__system_assembly_base_h__
#define __mast__system_assembly_base_h__

// C++ includes
#include <map>
#include <memory>


// MAST includes
#include "base/mast_data_types.h"


// libMesh includes
#include "libmesh/system.h"

namespace MAST {
    
    // Forward declerations
    class PhysicsDisciplineBase;
    class SystemInitialization;
    class ElementBase;
    class OutputFunctionBase;
    class MeshFieldFunction;
    class NonlinearSystem;
    class FEBase;
    class AssemblyElemOperations;
    
    class AssemblyBase {
    public:
        
        /*!
         *  constructor takes a reference to the discipline that provides
         *  the boundary conditions, volume loads, properties, etc.
         */
        AssemblyBase();
        
        /*!
         *   virtual destructor
         */
        virtual ~AssemblyBase();
        
        
        class SolverMonitor {
        public:
            SolverMonitor(){}
            virtual ~SolverMonitor(){}
            virtual void init(MAST::AssemblyBase& assembly) = 0;
            virtual void clear() = 0;
        };
        
        /*!
         *   @returns a const reference to the PhysicsDisciplineBase object
         *   associated with this object
         */
        const MAST::PhysicsDisciplineBase& discipline() const;
        
        /*!
         *   @returns a non-const reference to the PhysicsDisciplineBase object
         *   associated with this object
         */
        MAST::PhysicsDisciplineBase& discipline();
        
        
        /*!
         *   @returns a reference to the element operations object
         */
        MAST::AssemblyElemOperations& get_elem_ops();
        
        
        /*!
         *   attaches a system to this discipline, and vice-a-versa
         */
        virtual void
        attach_discipline_and_system(MAST::AssemblyElemOperations& elem_ops,
                                     MAST::PhysicsDisciplineBase& discipline,
                                     MAST::SystemInitialization& system);

        
        /*!
         *   Reattaches to the same system that was attached earlier.
         *
         *   This cannot be called if the clear_discipline_and_system() method
         *   has been called.
         */
        virtual void
        reattach_to_system() = 0;

        
        /*!
         *   clears association with a system to this discipline, and vice-a-versa
         */
        virtual void
        clear_discipline_and_system( );

        
        /*!
         *   @returns a const reference to the MAST::NonlinearSystem object
         *   associated with this object
         */
        const MAST::NonlinearSystem& system() const;
        
        /*!
         *   @returns a non-const reference to the MAST::NonlinearSystem object
         *   associated with this object
         */
        MAST::NonlinearSystem& system();

        /*!
         *   @returns a non-const reference to the MAST::SystemInitialization
         */
        MAST::SystemInitialization& system_init();


        /*!
         *   attaches the solver monitor, which is a user provided routine
         *   that is called each time
         */
        void set_solver_monitor(MAST::AssemblyBase::SolverMonitor& monitor);

        /*!
         *   @returns the solver monitor, which is a user provided routine
         *   that is called each time
         */
        MAST::AssemblyBase::SolverMonitor* get_solver_monitor();

        /*!
         *   clears the monitor object
         */
        void clear_solver_monitor();
        
        /*!
         *   tells the assembly object that this function is will
         *   need to be initialized before each residual evaluation
         */
        void attach_solution_function(MAST::MeshFieldFunction& f);

        
        /*!
         *   removes the attachment of the solution function
         */
        void detach_solution_function();
        
        
        /*!
         *   evaluates the volume and boundary outputs for the specified
         *   solution
         */
        virtual void calculate_outputs(const libMesh::NumericVector<Real>& X);

        
        /*!
         *   evaluates the sensitivity of the outputs in the attached 
         *   discipline with respect to the parametrs in \par params.
         *   The base solution should be provided in \par X. If the parameter 
         *   \par if_total_sensitivity is true, then the method will calculate 
         *   the total derivative of the output with respect to the parameters.
         *   In this case, the method will look for the sensitivity solution for 
         *   the i^th parameter in get_sensitivity_solution() of the associated
         *   libMesh::System object. If the parameter \par if_total_sensitivity
         *   if \p false, then the method will calculate the sensitivity
         *   assuming the sensitivity of solution is zero. This can be used for
         *   adjoint sensitivity analysis.
         */
        void calculate_output_sensitivity(libMesh::ParameterVector& params,
                                          const bool if_total_sensitivity,
                                          const libMesh::NumericVector<Real>& X);

        /*!
         *   localizes the parallel vector so that the local copy
         *   stores all values necessary for calculation of the
         *   element quantities
         */
        std::unique_ptr<libMesh::NumericVector<Real> >
        build_localized_vector(const libMesh::System& sys,
                               const libMesh::NumericVector<Real>& global);
        
        
        /*!
         *   @returns a MAST::FEBase object for calculation of finite element
         *   quantities. For all standard applications this is a wrapper
         *   around the libMesh::FEBase class, which is specialized for
         *   cut-cell applications where a sub-finite element is created
         *   for element integration.
         */
        virtual std::unique_ptr<MAST::FEBase>
        build_fe(const libMesh::Elem& e);
        

    protected:
        
        
        
        
        /*!
         *   assembles the outputs for this element
         */
        virtual void
        _elem_outputs(MAST::ElementBase& elem,
                      std::multimap<libMesh::subdomain_id_type, MAST::OutputFunctionBase *> &vol_output,
                      std::multimap<libMesh::boundary_id_type,MAST::OutputFunctionBase *> &side_output);

        
        /*!
         *   assembles the sensitivity of outputs for this element
         */
        virtual void
        _elem_output_sensitivity(MAST::ElementBase& elem,
                                 std::multimap<libMesh::subdomain_id_type, MAST::OutputFunctionBase *> &vol_output,
                                 std::multimap<libMesh::boundary_id_type,MAST::OutputFunctionBase *> &side_output);

        /*!
         *   provides assembly elem operations for use by this class
         */
        MAST::AssemblyElemOperations*  _elem_ops;

        
        /*!
         *   PhysicsDisciplineBase object for which this class is assembling
         */
        MAST::PhysicsDisciplineBase* _discipline;
        
        
        /*!
         *   System for which this assembly is performed
         */
        MAST::SystemInitialization* _system;
        
        /*!
         *   system solution that will be initialized before each solution
         */
        MAST::MeshFieldFunction* _sol_function;
        
        /*!
         *   User provided solver monitor is attached to the linear
         *   nonlinear solvers, if provided
         */
        MAST::AssemblyBase::SolverMonitor *_solver_monitor;
    };
        
}

#endif // __mast__system_assembly_base_h__

