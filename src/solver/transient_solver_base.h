

#ifndef __mast__transient_solver_base__
#define __mast__transient_solver_base__

// MAST includes
#include "base/mast_data_types.h"

// libMesh includes
#include "libmesh/transient_system.h"
#include "libmesh/nonlinear_implicit_system.h"


namespace MAST {
    
    // Forward declerations
    class TransientAssembly;
    class ElementBase;
    
    class TransientSolverBase {
    public:
        TransientSolverBase();
        
        virtual ~TransientSolverBase();

        /*!
         *   Attaches the assembly object that provides the x_dot, M and J
         *   quantities for the element
         */
        void set_assembly(MAST::TransientAssembly& assembly);

        /*!
         *   Clears the assembly object
         */
        void clear_assembly();
        
        /*!
         *   time step
         */
        Real dt;

        /*!
         *    @returns a const reference to the localized solution from 
         *    iteration number = current - prev_iter. So, \par prev_iter = 0
         *    gives the current solution estimate. Note that \par prev_iter
         *    cannot be greater than the total number of iterations that this
         *    solver stores solutions for.
         */
        const libMesh::NumericVector<Real>&
        solution(unsigned int prev_iter = 0) const;
        
        /*!
         *    @returns a const reference to the localized velocity from
         *    iteration number = current - prev_iter. So, \par prev_iter = 0
         *    gives the current velocity estimate. Note that \par prev_iter
         *    cannot be greater than the total number of iterations that this
         *    solver stores solutions for.
         */
        const libMesh::NumericVector<Real>&
        velocity(unsigned int prev_iter = 0) const;

        
        /*!
         *   solves the current time step for solution and velocity
         */
        virtual void solve() = 0;
        
        
        /*!
         *   advances the time step and copies the current solution to old
         *   solution, and so on.
         */
        virtual void advance_time_step();

        
        /*!
         *    TransientAssembly needs to be able to call the assembly routines
         *    of this class.
         */
        friend class MAST::TransientAssembly;

        
    protected:
        
        /*!
         *    @returns the number of iterations for which solution and velocity
         *    are to be stored.
         */
        virtual unsigned int _n_iters_to_store() const = 0;

        /*!
         *    localizes the relevant solutions for system assembly.
         */
        virtual void _localize_solutions();
        
        /*!
         *   performs the element calculations over \par elem, and returns
         *   the element vector and matrix quantities in \par mat and
         *   \par vec, respectively. \par if_jac tells the method to also
         *   assemble the Jacobian, in addition to the residual vector.
         */
        virtual void
        _elem_calculations(MAST::ElementBase& elem,
                           const std::vector<libMesh::dof_id_type>& dof_indices,
                           bool if_jac,
                           RealVectorX& vec,
                           RealMatrixX& mat) = 0;
        
        /*!
         *   performs the element sensitivity calculations over \par elem,
         *   and returns the element residual sensitivity in \par vec .
         */
        virtual void
        _elem_sensitivity_calculations(MAST::ElementBase& elem,
                                       const std::vector<libMesh::dof_id_type>& dof_indices,
                                       RealVectorX& vec) = 0;
        
        /*!
         *   Associated TransientAssembly object that provides the 
         *   element level quantities
         */
        MAST::TransientAssembly* _assembly;
        
        /*!
         *   TransientNonlinearImplicitSystem for which this object is
         *   calculating the solution
         */
        libMesh::TransientNonlinearImplicitSystem* _system;
        
        /*!
         *   localized solution vector from previous time step needed for
         *   element calculations on local processor
         */
        std::vector<libMesh::NumericVector<Real>*> _solution;

        
        /*!
         *   localized solution vector needed for element calculations
         *   on local processor
         */
        std::vector<libMesh::NumericVector<Real>*> _velocity;
        
    };

}


#endif // __mast__transient_solver_base__
