#ifndef __mast__elem_base__
#define __mast__elem_base__

// MAST includes
#include "base/mast_data_types.h"

// libMesh includes
#include "libmesh/elem.h"
#include "libmesh/system.h"
#include "libmesh/fe_base.h"
#include "libmesh/quadrature.h"



namespace MAST {
    
    // Forward declerations
    class FunctionBase;
    class SystemInitialization;
    class LocalElemBase;
    
    
    /*!
     *    This is the base class for elements that implement calculation of
     *    finite element quantities over the domain and sides of a geometric
     *    element.
     *
     *    The primary functions required of the elements depend on the nature
     *    of analysis:
     *    - Nonlinear Analysis
     *    \f[ r(X,p) = 0 \f]
     *    A nonliner analysis requires element residuals and Jacobians.
     *    Sensitivity analysis of the nonlinear problem requires
     *    the sensitivity of the residual vector with respect to the concerned
     *    parameter.
     *
     *    - Eigenvalue Analysis
     *    \f[ \lambda A(X,p) = B(X,p) \f]
     *    An eigenvalue analysis requires the coefficent matrices \f$A(X,p)\f$
     *    and \f$B(X,p)\f$. For cases where the eigenvalue problem is defined
     *    using small disturbances about a steady-state solution, the base
     *    solution (and its sensitivity for sensitivity problems) needs to be
     *    provided to the element.
     *
     *    - Transient Analysis for first order systems
     *    \f[ m(\dot{X}, X, p) = f(X,p) \f]
     *
     */
    class ElementBase {
    public:
        
        /*!
         *   The default constructor.
         *   \param sys SystemInitialization object which provides the system
         *   for which this element will perorm the calculations.
         *   \param elem libMesh::Elem object on which calculations will be
         *   performed.
         */
        ElementBase(MAST::SystemInitialization& sys,
                    const libMesh::Elem& elem);
        
        
        /*!
         *   Default virtual destructor
         */
        virtual ~ElementBase();
        
        
        /*!
         *   @returns a reference to the libMesh::System object
         */
        libMesh::System& system();
        
        
        /*!
         *   @returns a constant reference to the element
         */
        const libMesh::Elem& elem() const {
            return _elem;
        }
        
        
        /*!
         *   @returns a constant reference to the element in the local
         *   coordinate. This is needed for 1D or 2D elements that live
         *   in a 3D space.
         */
        const MAST::LocalElemBase& local_elem() const {
            return *_local_elem;
        }
        
        
        
        /*!
         *   @returns a constant reference to the volume quadrature rule
         */
        libMesh::QBase& quadrature_rule()  {
            return *_qrule;
        }
        
        
        
        /*!
         *   @returns a constant reference to the finite element object
         */
        libMesh::FEBase& fe()  {
            return *_fe;
        }
        
        
        /*!
         *   @returns a const reference to the solution vector, or its
         *   sensitivity if \p if_sens is true.
         */
        const RealVectorX& sol(bool if_sens = false) const;
        
        
        /*!
         *   stores \p vec as solution for element level calculations,
         *   or its sensitivity if \p if_sens is true.
         */
        virtual void set_solution(const RealVectorX& vec,
                                  bool if_sens = false);
        
        
        /*!
         *    stores \p vec as velocity for element level calculations,
         *    or its sensitivity if \p if_sens is true.
         */
        virtual void set_velocity(const RealVectorX& vec,
                                  bool if_sens = false);

        
        /*!
         *    stores \p vec as acceleration for element level calculations,
         *    or its sensitivity if \p if_sens is true.
         */
        virtual void set_acceleration(const RealVectorX& vec,
                                      bool if_sens = false);

        
        /*!
         *   This is used for cases where a linearized problem is solved
         *   about a stationary base solution. This method stores
         *   \p vec as the base solution, or its sensitivity if \p
         *   if_sens is true.
         */
        virtual void set_base_solution(const RealVectorX& vec,
                                       bool if_sens = false);
        
        /*!
         *   Attaches the function that represents the system solution
         */
        void attach_active_solution_function(MAST::FunctionBase &f);
        
        
        /*!
         *   Detaches the function object that may have been attached to the
         *   element.
         */
        void detach_active_solution_function();

        
        /*!
         *   parameter for which sensitivity has to be calculated.
         */
        const MAST::FunctionBase* sensitivity_param;
        
        
        /*!
         *   @returns a constant reference to the geometric element used for
         *   initialization of finite element quadrature and shape functions.
         *   This is needed for cases where a 1D or 2D element might live in a
         *   3D space, in which case the element returned will be one that
         *   has been transformed to a local coordinate system. For a 3D element,
         *   the method returns the element used to initialize this object.
         */
        const libMesh::Elem& get_elem_for_quadrature() const;
        
    protected:
        
        
        /*!
         *   Initializes the quadrature and finite element for element volume
         *   integration.
         *   \param e libMesh::Elem for which the finite element is initialized
         */
        virtual void
        _init_fe_and_qrule(const libMesh::Elem& e);
        
        
        /*!
         *   @returns the quadrature and finite element for element side
         *   integration. These are raw pointers created using new. The
         *   pointers must be deleted at the end of scope.
         */
        virtual void
        _get_side_fe_and_qrule(const libMesh::Elem& e,
                               unsigned int s,
                               std::auto_ptr<libMesh::FEBase>& fe,
                               std::auto_ptr<libMesh::QBase>& qrule);
        
        
        /*!
         *   SystemInitialization object associated with this element
         */
        MAST::SystemInitialization& _system;
        
        
        /*!
         *   geometric element for which the computations are performed
         */
        const libMesh::Elem& _elem;
        
        
        /*!
         *   pointer to the active solution mesh field function. If this 
         *   has been set, then some of the element properties are
         *   dependent on the element solution, and the element should 
         *   perform the necessary operations in calculation of the Jacobian
         */
        MAST::FunctionBase* _active_sol_function;
        
        
        /*!
         *   local element to support the presence of 1D and 2D elements
         *   in 3D space
         */
        std::auto_ptr<MAST::LocalElemBase> _local_elem;
        
        
        /*!
         *    time for which system is being assembled
         */
        const Real& _time;
        
        
        /*!
         *   local solution
         */
        RealVectorX _sol;
        
        
        /*!
         *   local solution sensitivity
         */
        RealVectorX _sol_sens;
        
        
        /*!
         *   local velocity
         */
        RealVectorX _vel;
        
        
        /*!
         *   local velocity
         */
        RealVectorX _vel_sens;

        
        /*!
         *   local acceleration
         */
        RealVectorX _accel;
        
        
        /*!
         *   local acceleration
         */
        RealVectorX _accel_sens;
        
        
        /*!
         *   base solution about which a linearized solution is performed
         */
        RealVectorX _base_sol;
        
        
        /*!
         *   base solution sensitivity
         */
        RealVectorX _base_sol_sens;
        
        
        /*!
         *   element finite element for computations
         */
        std::auto_ptr<libMesh::FEBase> _fe;
        
        
        /*!
         *   element quadrature rule for computations
         */
        std::auto_ptr<libMesh::QBase> _qrule;
    };
}


#endif // __mast__elem_base__
