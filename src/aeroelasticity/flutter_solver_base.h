/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2017  Manav Bhatia
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

#ifndef __mast__flutter_solver_base_h__
#define __mast__flutter_solver_base_h__

// C++ includes
#include <string>
#include <fstream>
#include <iomanip>


// MAST includes
#include "base/mast_data_types.h"


// libMesh includes
#include "libmesh/parameter_vector.h"
#include "libmesh/numeric_vector.h"



namespace MAST {
    
    // Forward declerations
    class Parameter;
    class FlutterModel;
    class FlutterRootBase;
    class FlutterSolutionBase;
    class FlutterRootCrossoverBase;
    class StructuralFluidInteractionAssembly;
    template <typename ValType> class BasisMatrix;
    
    
    class FlutterSolverBase {
        
    public:
        
        /*!
         *   defalut constructor
         */
        FlutterSolverBase();
        
        
        virtual ~FlutterSolverBase();
        
        
        /*!
         *   abstract class defines the interface to provide the steady-state
         *   solution
         */
        class SteadySolver {
        public:
            
            
            SteadySolver() {}
            
            virtual ~SteadySolver() {}
            
            /*!
             *  solves for the steady state solution, and @returns
             *  a const-reference to the solution.
             */
            virtual const libMesh::NumericVector<Real>&
            solve() = 0;

            
            /*!
             * @returns  a const-reference to the solution.
             */
            virtual const libMesh::NumericVector<Real>&
            solution() const = 0;
            
        };

        
        /*!
         *    attaches the assembly object to this solver.
         */
        void attach_assembly(MAST::StructuralFluidInteractionAssembly&   assembly);
        
        
        /*!
         *    attaches the steady solution object
         */
        void attach_steady_solver(MAST::FlutterSolverBase::SteadySolver& solver);
        
        
        /*!
         *   clears the solution and other data from this solver
         */
        virtual void clear();
        
        
        /*!
         *   clears the assembly object
         */
        virtual void clear_assembly_object();
        
        
        /*!
         *    initializes the data structres for a flutter solution.
         */
        void initialize(std::vector<libMesh::NumericVector<Real>*>& basis);

        
        
        void set_output_file(const std::string& nm) {
            
            if (!_output)
                _output = new std::ofstream;
            
            _output->close();
            _output->open(nm.c_str(), std::ofstream::out);
        }
        
        
        /*!
         *   Prints the sorted roots to the \par output
         */
        virtual void print_sorted_roots() = 0;
        
        
        /*!
         *   Prints the crossover points output. If no pointer to output is given
         *   then the output defined by set_output_file() is used.
         */
        virtual void print_crossover_points() = 0;
        
        
    protected:
        
        
        /*!
         *   structural assembly that provides the assembly of the system
         *   matrices.
         */
        MAST::StructuralFluidInteractionAssembly*       _assembly;
        
        
        /*!
         *   basis vector used to define the reduced order model
         */
        std::vector<libMesh::NumericVector<Real>*>*     _basis_vectors;
        
        
        /*!
         *    file to which the result will be written
         */
        std::ofstream*                                  _output;
        
        
        /*!
         *    object provides the steady state solution.
         */
        MAST::FlutterSolverBase::SteadySolver* _steady_solver;
        
    };
}


#endif // __mast__flutter_solver_base_h__
