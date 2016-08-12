/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2016  Manav Bhatia
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

#ifndef __mast__multiphysics_nonlinear_solver_base_h__
#define __mast__multiphysics_nonlinear_solver_base_h__

// C++ includes
#include <vector>
#include <string>

// PETSc includes
#include <petscmat.h>


namespace MAST {

    // Forward declerations
    class NonlinearImplicitAssembly;
    
    
    class MultiphysicsNonlinearSolverBase {
        
    public:
        
        /*!
         *   default constructor
         */
        MultiphysicsNonlinearSolverBase(const std::string& nm,
                                        unsigned int n);
        

        /*!
         *   destructor
         */
        virtual ~MultiphysicsNonlinearSolverBase();

        
        /*!
         *    @returns the name of this solver
         */
        const std::string name() const {
            
            return _name;
        }
        
        
        /*!
         *   @returns the number of systems
         */
        unsigned int n_disciplines() const {
            
            return _n_disciplines;
        }
        
        
        
        /*!
         *   method to set the n^th discipline of this multiphysics system
         *   assembly.
         */
        void set_system_assembly(unsigned int i,
                                 MAST::NonlinearImplicitAssembly& assembly);


        
        /*!
         *   @returns a reference to the n^th discipline of this multiphysics
         *   system assembly
         */
        MAST::NonlinearImplicitAssembly& get_system_assembly(unsigned int i);
        

        /*!
         *   @returns the MPI communicator context for the global system
         */
        MPI_Comm comm() {
            
            return _g_comm;
        }
        

        /*!
         *   @returns a reference to the petsc index sets
         */
        std::vector<IS>& index_sets() {
            
            return _is;
        }
        
        
        /*!
         *   @returns the global Jacobian matrix context
         */
        Mat mat() {
            
            return _mat;
        }
        
        /*!
         *   solves the system using the nested matrices that uses the
         *   discipline specific solver options
         */
        void solve();
        

        
    protected:
        
        
        /*!
         *  name of this multiphysics solution
         */
        const std::string  _name;
        
        
        /*!
         *   number of disciplines
         */
        const unsigned int _n_disciplines;
        
        
        /*!
         *   vector of assembly objects for each discipline in this 
         *   multiphysics system
         */
        std::vector<MAST::NonlinearImplicitAssembly*>  _discipline_assembly;

        MPI_Group        _g_union;
        MPI_Comm         _g_comm;
        
        std::vector<IS>  _is;
        std::vector<Mat> _sub_mats; // row-major ordering
        unsigned int     _n_dofs;

        Mat              _mat;
        Vec              _sol, _res;

    };
}



#endif // __mast__multiphysics_nonlinear_solver_base_h__
