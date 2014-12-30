
#ifndef __mast__first_order_newmark_transient_solver__
#define __mast__first_order_newmark_transient_solver__

// MAST includes
#include "solver/transient_solver_base.h"


namespace MAST {
    
    
    /*!
     *    This class implements the Newmark solver for solution of a
     *    first-order ODE.
     *
     *    The system is represented as 
     *    \f[ x^k_t = x_{t_0} + (1-\beta)dt\dot{x}_{t_0} + \beta dt\dot{x}^k_t \f],
     *    where \f$k\f$ is the nonlinear iteration for each time step from 
     *    \f$t_0\f$ to \f$t\f$, and $\beta=[0,1]$ is a user defined factor.
     *
     *    A finite element solver provides the following semi-discrete form for
     *    a first order transient system:
     *    \f[ M^k \dot{x}^k_t + f(x^k_t) = 0 \f].
     *    Then, to solve for \f$x^k_t\f$, the expression for \f$x^k_t\f$ is 
     *    multiplied by \f$M^k\f$ and written in a residual form as
     *    \f[ r(x^k_t) = M^k \left( x^k_t - x_{t_0} - (1-\beta)dt\dot{x}_{t_0}\right) + \beta dt f(x^k_t) \f]
     *    where the Jacobian is defined as 
     *    \f[ J^k = \frac{\partial r(x^k_t)}{\partial x^k_t} = M^k + \beta dt \frac{\partial f(x^k_t)}{\partial x^k_t} \f].
     *    This assumes that the mass matrix is independent of the \f$x^k_t\f$.
     *
     */
    class FirstOrderNewmarkTransientSolver:
    public MAST::TransientSolverBase {
    public:
        FirstOrderNewmarkTransientSolver();
        
        virtual ~FirstOrderNewmarkTransientSolver();
        
        /*!
         *    \f$ \beta \f$ parameter used by this solver.
         */
        Real beta;
        
        /*!
         *   solves the current time step for solution and velocity
         */
        virtual void solve();
        
        /*!
         *   advances the time step and copies the current solution to old
         *   solution, and so on.
         */
        virtual void advance_time_step();
        
    protected:
        
        /*!
         *    @returns the number of iterations for which solution and velocity
         *    are to be stored.
         */
        virtual unsigned int _n_iters_to_store() const {
            return 2;
        }
        
        /*!
         *    provides the element with the transient data for calculations
         */
        virtual void _set_element_data(std::vector<libMesh::dof_id_type>& dof_indices,
                                       MAST::ElementBase& elem);

        /*!
         *    update the transient solution based on the current solution
         */
        virtual void _update_velocity(libMesh::NumericVector<Real>& vec);
        
        
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
                           RealMatrixX& mat);
        
        /*!
         *   performs the element sensitivity calculations over \par elem,
         *   and returns the element residual sensitivity in \par vec .
         */
        virtual void
        _elem_sensitivity_calculations(MAST::ElementBase& elem,
                                       const std::vector<libMesh::dof_id_type>& dof_indices,
                                       RealVectorX& vec);
    };
    
}

#endif // __mast__first_order_newmark_transient_solver__
