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

#ifndef __mast__elem_base__
#define __mast__elem_base__

// C++ includes
#include <map>

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
    class OutputFunctionBase;
    
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
        
        /*!
         *   evaluates an output quantity requested in the map over the
         *   boundary of the element that may coincide with the boundary
         *   identified in the map. The derivative with respect to the
         *   state variables is provided if \p request_derivative is true.
         */
        virtual bool
        volume_output_quantity (bool request_derivative,
                                bool request_sensitivity,
                                std::multimap<libMesh::subdomain_id_type, MAST::OutputFunctionBase*>& output) = 0;
        
        
        /*!
         *   evaluates an output quantity requested in the map over the
         *   boundary of the element that may coincide with the boundary
         *   identified in the map. The derivative with respect to the
         *   state variables is provided if \p request_derivative is true.
         */
        virtual bool
        side_output_quantity (bool request_derivative,
                              std::multimap<libMesh::boundary_id_type, MAST::OutputFunctionBase*>& output) = 0;

    
    protected:
        
        
        /*!
         *   Initializes the quadrature and finite element for element volume
         *   integration.
         *   \param e libMesh::Elem for which the finite element is initialized.
         *   \param pts the points at which the element should be initialized. If NULL, 
         *    the points specified by quadrature rule will be used.
         */
        virtual void
        _init_fe_and_qrule(const libMesh::Elem& e,
                           libMesh::FEBase **fe,
                           libMesh::QBase **qrule,
                           const std::vector<libMesh::Point>* pts = NULL);
        
        
        /*!
         *   @returns the quadrature and finite element for element side
         *   integration. These are raw pointers created using new. The
         *   pointers must be deleted at the end of scope. The last argument
         *   \p if_calculate_dphi tell the function to request the 
         *   \p fe object to also initialize the calculation of shape function
         *   derivatives 
         */
        virtual void
        _get_side_fe_and_qrule(const libMesh::Elem& e,
                               unsigned int s,
                               libMesh::FEBase **fe,
                               libMesh::QBase **qrule,
                               bool if_calculate_dphi);
        
        
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
        libMesh::FEBase *_fe;
        
        
        /*!
         *   element quadrature rule for computations
         */
        libMesh::QBase *_qrule;
    };
}


#endif // __mast__elem_base__
