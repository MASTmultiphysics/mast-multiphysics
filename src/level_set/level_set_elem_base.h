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


#ifndef __mast__level_set_elem_base__
#define __mast__level_set_elem_base__

// MAST includes
#include "base/elem_base.h"

namespace MAST {
    
    // Forward declerations
    class ElementPropertyCardBase;
    class BoundaryConditionBase;
    class FEMOperatorMatrix;
    template <typename ValType> class FieldFunction;
    
    
    class LevelSetElementBase:
    public MAST::ElementBase
    {
    public:
        /*!
         *   Constructor
         */
        LevelSetElementBase(MAST::SystemInitialization&             sys,
                            MAST::AssemblyBase&                     assembly,
                            const MAST::GeomElem&                   elem);
        
        
        virtual ~LevelSetElementBase();
        
        
        void
        set_velocity_function(const MAST::FieldFunction<Real>& vel) {
            
            libmesh_assert(!_phi_vel);
            _phi_vel = &vel;
        }
        
        /*!
         *   This can operate in one of two modes: propagation of level set
         *   given Vn, or reinitialization of level set so that |grad(phi)|=1.
         *   This method sets the flag for propagation to \p true or \p false.
         */
        void set_propagation_mode(bool f) {
            _if_propagation = f;
        }

        
        /*!
         *   For reinitialization to \f$ |\nabla(\phi)| = 1 \f$, the solution
         *   before initialization is used to calculate the source and velocity
         *   switching. This method sets that solution
         */
        void set_reference_solution_for_initialization(const RealVectorX& sol);
        
        
        /*!
         *   internal force contribution to system residual
         */
        virtual bool
        internal_residual (bool request_jacobian,
                           RealVectorX& f,
                           RealMatrixX& jac);
        
        
        /*!
         *   inertial force contribution to system residual
         */
        virtual bool
        velocity_residual (bool request_jacobian,
                           RealVectorX& f,
                           RealMatrixX& jac_xdot,
                           RealMatrixX& jac);
        
        /*!
         *   side external force contribution to system residual
         */
        bool
        side_external_residual (bool request_jacobian,
                                RealVectorX& f,
                                RealMatrixX& jac,
                                std::multimap<libMesh::boundary_id_type, MAST::BoundaryConditionBase*>& bc);
        
        /*!
         *   volume external force contribution to system residual
         */
        bool
        volume_external_residual (bool request_jacobian,
                                  RealVectorX& f,
                                  RealMatrixX& jac,
                                  std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*>& bc);
        
        /*!
         *   sensitivity of the internal force contribution to system residual
         */
        virtual bool
        internal_residual_sensitivity (bool request_jacobian,
                                       RealVectorX& f,
                                       RealMatrixX& jac);
        /*!
         *   sensitivity of the damping force contribution to system residual
         */
        virtual bool
        velocity_residual_sensitivity (bool request_jacobian,
                                       RealVectorX& f,
                                       RealMatrixX& jac);
        
        /*!
         *   sensitivity of the side external force contribution to system residual
         */
        bool
        side_external_residual_sensitivity (bool request_jacobian,
                                            RealVectorX& f,
                                            RealMatrixX& jac,
                                            std::multimap<libMesh::boundary_id_type, MAST::BoundaryConditionBase*>& bc);
        
        /*!
         *   sensitivity of the volume external force contribution to system residual
         */
        bool
        volume_external_residual_sensitivity (bool request_jacobian,
                                              RealVectorX& f,
                                              RealMatrixX& jac,
                                              std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*>& bc);
        
        /*!
         *   @returns the element volume over the domain of the element that has
         *   positive level set function
         */
        Real
        volume();

        /*!
         *   Approximates the integral of the Dirac delta function to approximate
         *   the perimeter. The approximation of Dirac delta function is obtained
         *   from the derivative of an approximation to the Heaviside function
         *   \f[ H_d(\phi) \approx \frac{1}{2} \left( 1 + \frac{2}{\pi} \arctan(\phi/d) \right).  \f]
         *   Then, the derivative defines the Dirac delta function
         *   \f[ \delta_d(\phi) = \frac{dH_d(\phi)}{d\phi} = \frac{1}{\pi d} \frac{1}{1+(\phi/d)^2} .  \f]
         *   Sensitivity analysis requires the derivative of this function with
         *   respect to a variable, which is expressed as
         *   \f[ \frac{d\delta_d(\phi)}{d\alpha} = \frac{-2 \phi}{\pi d^3} \frac{1}{\left(1+(\phi/d)^2\right)^2} \frac{d\phi}{d\alpha}.  \f]
         *
         *   @returns the computed integral of \f$ \delta_d(\phi) \f$
         */
        Real
        perimeter();

        /*!
         *   computes the partial derivative of the integral of the Dirac
         *   delta function using the solution and sensitivity solution
         *   set for this element. 
         */
        Real
        perimeter_sensitivity();

        
        /*!
         *   @returns the contribution of the side to
         *   \f$ \int_\Gamma V_n d\Gamma \f$, where \f$ V_n \f$ is the boundary
         *   normal velocity evaluated from the solution and sensitivity.
         */
        Real
        volume_boundary_velocity_on_side(unsigned int s);

    protected:

        /*!
         *   calculates the velocity at the quadrature point
         */
        void _velocity_and_source(const unsigned int qp,
                                  const libMesh::Point& p,
                                  const Real t,
                                  const MAST::FEMOperatorMatrix& Bmat,
                                  const std::vector<MAST::FEMOperatorMatrix>& dBmat,
                                  RealVectorX& vel,
                                  Real&        source);
        
        /*!
         *   initializes the tau operator
         */
        void _tau(const MAST::FEBase& fe,
                  unsigned int qp,
                  const MAST::FEMOperatorMatrix& Bmat,
                  const std::vector<MAST::FEMOperatorMatrix>& dBmat,
                  const RealVectorX& vel,
                  RealMatrixX& tau);
        
        
        
        
        void
        _dc_operator(const MAST::FEBase& fe,
                     const unsigned int qp,
                     const std::vector<MAST::FEMOperatorMatrix>& dB_mat,
                     const RealVectorX& vel,
                     Real& dc);
        
        void
        _calculate_dxidX (const MAST::FEBase& fe,
                          const unsigned int qp,
                          RealMatrixX& dxi_dX);
        
        /*!
         *    When \p mass = false, initializes the FEM operator matrix to the
         *    shape functions as
         *    \f[  B = \left[ \begin{array}{c}
         *    {\bf N} \\ {\bf N} \\ {\bf N}
         *    \end{array} \right] \f]
         *    \f[  dB[0] = \frac{\partial {\bf N}}{\partial x} \f]
         *
         *    \f[  dB[1] = \frac{\partial {\bf N}}{\partial y} \f]
         *
         *    \f[  dB[2] = \frac{\partial {\bf N}}{\partial z} \f]
         */
        void _initialize_fem_operators(const unsigned int qp,
                                       const MAST::FEBase& fe,
                                       MAST::FEMOperatorMatrix& Bmat,
                                       std::vector<MAST::FEMOperatorMatrix>& dBmat);
        
        
        /*!
         *   element property
         */
        const MAST::FieldFunction<Real>*  _phi_vel;

        /*!
         *   this can operate in one of two modes: propagation of level set
         *   given Vn, or reinitialization of level set so that |grad(phi)|=1
         */
        bool  _if_propagation;
        
        /*!
         *     reference solution for reinitialization of the level set
         */
        RealVectorX    _ref_sol;
    };
    
}


#endif // __mast__level_set_elem_base__

