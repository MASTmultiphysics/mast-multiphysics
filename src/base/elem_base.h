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

#ifndef __mast__elem_base__
#define __mast__elem_base__

// C++ includes
#include <map>
#include <memory>

// MAST includes
#include "base/mast_data_types.h"

// libMesh includes
#include "libmesh/elem.h"
#include "libmesh/system.h"



namespace MAST {
    
    // Forward declerations
    class FunctionBase;
    class SystemInitialization;
    class GeomElem;
    class NonlinearSystem;
    class FEBase;
    class AssemblyBase;
    
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
     *    \f[ A(X,p) = \lambda  B(X,p) \f]
     *    An eigenvalue analysis requires the coefficent matrices \f$A(X,p)\f$
     *    and \f$B(X,p)\f$. For cases where the eigenvalue problem is defined
     *    using small disturbances about a steady-state solution, the base
     *    solution (and its sensitivity for sensitivity problems) needs to be
     *    provided to the element.
     *
     *    - Transient Analysis for first order systems
     *    \f[ f_m(\dot{X}, X, p) + f(X,p) = 0 \f]
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
        ElementBase(MAST::SystemInitialization&     sys,
                    MAST::AssemblyBase&             assembly,
                    const MAST::GeomElem&           elem);
        
        
        /*!
         *   Default virtual destructor
         */
        virtual ~ElementBase();
        
        
        /*!
         *   @returns a reference to the libMesh::System object
         */
        MAST::SystemInitialization& system_initialization() {
            return _system;
        }

        /*!
         *   @returns a reference to the libMesh::System object
         */
        MAST::AssemblyBase& assembly() {
            return _assembly;
        }

        
        /*!
         *   @returns a reference to the libMesh::System object
         */
        MAST::NonlinearSystem& system();
        
        
        /*!
         *   @returns a constant reference to the element
         */
        const MAST::GeomElem& elem() const {
            return _elem;
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
         *   This provides the perturbed solution (or its sensitivity if
         *   \p if_sens is true.) for linearized analysis.
         *   stores \p vec as solution for element level calculations.
         */
        virtual void set_perturbed_solution(const RealVectorX& vec,
                                            bool if_sens = false);

        
        /*!
         *   This provides the complex solution (or its sensitivity if 
         *   \p if_sens is true.) for frequecy-domain analysis.
         *   stores \p vec as solution for element level calculations.
         */
        virtual void set_complex_solution(const ComplexVectorX& vec,
                                          bool if_sens = false);

        
        /*!
         *    stores \p vec as velocity for element level calculations,
         *    or its sensitivity if \p if_sens is true.
         */
        virtual void set_velocity(const RealVectorX& vec,
                                  bool if_sens = false);

        
        /*!
         *    stores \p vec as perturbed velocity for element level
         *    calculations, or its sensitivity if \p if_sens is true.
         */
        virtual void set_perturbed_velocity(const RealVectorX& vec,
                                            bool if_sens = false);

        
        
        /*!
         *    stores \p vec as acceleration for element level calculations,
         *    or its sensitivity if \p if_sens is true.
         */
        virtual void set_acceleration(const RealVectorX& vec,
                                      bool if_sens = false);

        
        /*!
         *    stores \p vec as perturbed acceleration for element level
         *    calculations, or its sensitivity if \p if_sens is true.
         */
        virtual void set_perturbed_acceleration(const RealVectorX& vec,
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
        
    
    protected:
        
        
        /*!
         *   SystemInitialization object associated with this element
         */
        MAST::SystemInitialization& _system;
        
        /*!
         *    Assembly object
         */
        MAST::AssemblyBase&        _assembly;
        
        /*!
         *   geometric element for which the computations are performed
         */
        const MAST::GeomElem&      _elem;
        
        /*!
         *   pointer to the active solution mesh field function. If this 
         *   has been set, then some of the element properties are
         *   dependent on the element solution, and the element should 
         *   perform the necessary operations in calculation of the Jacobian
         */
        MAST::FunctionBase* _active_sol_function;
        
        
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
         *   local solution used for frequency domain analysis
         */
        ComplexVectorX _complex_sol;

        
        /*!
         *   local solution used for frequency domain analysis
         */
        ComplexVectorX _complex_sol_sens;

        
        /*!
         *   local solution used for linearized analysis
         */
        RealVectorX _delta_sol;
        
        
        /*!
         *   local solution used for linearized analysis
         */
        RealVectorX _delta_sol_sens;

        /*!
         *   local velocity
         */
        RealVectorX _vel;
        
        
        /*!
         *   local velocity
         */
        RealVectorX _vel_sens;

        
        /*!
         *   local velocity
         */
        RealVectorX _delta_vel;
        
        
        /*!
         *   local velocity
         */
        RealVectorX _delta_vel_sens;

        
        /*!
         *   local acceleration
         */
        RealVectorX _accel;

        
        /*!
         *   local acceleration
         */
        RealVectorX _accel_sens;

        
        /*!
         *   local acceleration
         */
        RealVectorX _delta_accel;
        
        
        /*!
         *   local acceleration
         */
        RealVectorX _delta_accel_sens;

    };
}


#endif // __mast__elem_base__
