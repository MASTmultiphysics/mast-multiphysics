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

#ifndef __mast__stress_output_base__
#define __mast__stress_output_base__

// C++ includes
#include <map>
#include <vector>

// MAST includes
#include "base/mast_data_types.h"
#include "base/physics_discipline_base.h"
#include "base/output_assembly_elem_operations.h"


// libMesh includes
#include "libmesh/elem.h"

namespace MAST {


    // Forward declerations
    class FunctionBase;
    

    /*!
     *    Data structure provides the mechanism to store stress and strain
     *    output from a structural analysis. The user may specify one of three
     *    evaluation modes: centroid, element quadrature points, or user
     *    specified element local points.
     *
     *    This is, by default, a 3D tensor that should be initialized
     *    with 6 components. The first three components of stress are the direct
     *    stresses in the \f$ x, y, z\f$ directions, \f$ \sigma_{xx},
     *    \sigma_{yy}, \sigma_{zz}\f$, and the next three components are
     *     the shear stresses \f$ \tau_{xy}, \tau_{yz}, \tau_{xz} \f$.
     *
     *    The first three components of strain are the direct
     *    strains in the \f$x, y, z\f$ directions, \f$ \epsilon_{xx},
     *    \epsilon_{yy}, \epsilon_{zz} \f$, and the next three components
     *    are the engineering shear strains \f$ \gamma_{xy}, \gamma_{yz},
     *     \gamma_{xz} \f$.
     */
    class StressStrainOutputBase:
    public MAST::OutputAssemblyElemOperations {
        
    public:
    
        
        /*!
         *    This class provides a mechanism to store stress/strain values,
         *    their derivatives and sensitivity values corresponding to a 
         *    specific quadrature point on the element.
         */
        class Data {
            
        public:
            Data(const RealVectorX& stress,
                 const RealVectorX& strain,
                 const libMesh::Point& qp,
                 const libMesh::Point& xyz,
                 Real JxW);
 
            
            /*!
             *   @returns the point at which stress is evaluated, in the
             *   element coordinate system.
             */
            const libMesh::Point&
            point_location_in_element_coordinate() const;
            
            /*!
             *   @returns stress
             */
            const RealVectorX& stress() const;

            
            /*!
             *   @returns strain
             */
            const RealVectorX& strain() const;
            
            
            /*!
             *   @returns the quadrature point JxW
             */
            Real quadrature_point_JxW() const;
            
            /*!
             *   @returns von Mises stress
             */
            Real von_Mises_stress() const;

            
            /*!
             *   @returns derivative of von Mises stress wrt state vector
             */
            RealVectorX dvon_Mises_stress_dX() const;

            
            /*!
             *   @returns derivative of von Mises stress wrt sensitivity parameter
             */
            Real dvon_Mises_stress_dp(const MAST::FunctionBase* f) const;


            /*!
             *   adds the derivative data
             */
            void set_derivatives(const RealMatrixX& dstress_dX,
                                 const RealMatrixX& dstrain_dX);

            
            /*!
             *   @return the derivative data
             */
            const RealMatrixX& get_dstress_dX() const;

            
            /*!
             *   @return the derivative data
             */
            const RealMatrixX& get_dstrain_dX() const;

            
            /*!
             *   sets the sensitivity of the data with respect to a function
             */
            void set_sensitivity(const MAST::FunctionBase* f,
                                 const RealVectorX& dstress_df,
                                 const RealVectorX& dstrain_df);

            
            /*!
             *   @ returns the sensitivity of the data with respect to a 
             *   function
             */
            const RealVectorX&
            get_stress_sensitivity(const MAST::FunctionBase* f) const;

            
            /*!
             *   @ returns the sensitivity of the data with respect to a
             *   function
             */
            const RealVectorX&
            get_strain_sensitivity(const MAST::FunctionBase* f) const;

            
        protected:

            /*!
             *   stress data
             */
            RealVectorX                   _stress;
            
            
            /*!
             *  strain data
             */
            RealVectorX                   _strain;
            
            
            /*!
             *  map of sensitivity of the stress with respect to a parameter
             */
            std::map<const MAST::FunctionBase*, RealVectorX> _stress_sensitivity;

            
            /*!
             *  map of sensitivity of the strain with respect to a parameter
             */
            std::map<const MAST::FunctionBase*, RealVectorX> _strain_sensitivity;

            
            /*!
             *   derivative of stress wrt state vector
             */
            RealMatrixX   _dstress_dX;

            /*!
             *   derivative of strain data wrt state vector
             */
            RealMatrixX   _dstrain_dX;

            
            /*!
             *   quadrature point location in element coordinates
             */
            libMesh::Point            _qp;
            
            
            /*!
             *   quadrature point location in physical coordinates
             */
            libMesh::Point            _xyz;
            
            
            /*!
             *   quadrature point JxW (product of transformation Jacobian and
             *   quadrature weight) for use in definition of functionals
             */
            Real _JxW;
        };
        

        
        /*!
         *    default constructor
         */
        StressStrainOutputBase();
        
        virtual ~StressStrainOutputBase();
        
        
        /*!
         *   sets the \f$p-\f$norm for calculation of stress functional
         */
        void set_p_val(Real p) {
            
            _p_norm =  p;
        }
        
        /*!
         *   @returns the \f$p-\f$norm for calculation of stress functional
         */
        Real get_p_val() {
            
            return _p_norm;
        }
        
        /*!
         *   initialize for the element.
         */
        virtual void init(const libMesh::Elem& elem);

        /*!
         *   zeroes the output quantity values stored inside this object
         *   so that assembly process can begin.
         */
        virtual void zero();
        
        /*!
         *   @returns the output quantity value
         */
        virtual Real output();
        
        /*!
         *   @returns the output quantity sensitivity for parameter
         */
        virtual Real output_sensitivity(const MAST::FunctionBase& p);
        
        
        /*!
         *   @returns the output quantity derivative with respect to
         *   state vector
         */
        virtual RealVectorX output_derivative();
        
        
        /*!
         *    this is the abstract interface to be implemented by derived
         *    classes. This method calculates the quantity \f$ q(X, p) \f$.
         */
        virtual void evaluate();
        
        
        /*!
         *    this is the abstract interface to be implemented by derived
         *    classes. This method calculates the partial derivative of quantity
         *    \f[ \frac{\partial q(X, p)}{\partial p} \f]  with
         *    respect to parameter \f$ p \f$.
         */
        virtual void evaluate_sensitivity(const MAST::FunctionBase& p);
        
        
        /*!
         *    this is the abstract interface to be implemented by derived
         *    classes. This method calculates the quantity
         *    \f[ \frac{\partial q(X, p)}{\partial X} \f] for this
         *    output function.
         */
        virtual void evaluate_derivative();
                
        
        /*!
         *   some simulations frequently deal with 1D/2D elements in 3D space,
         *   which requires use of MAST::LocalElemFE.
         */
        virtual bool
        if_use_local_elem() const {
            
            return true;
        }
        
        
        /*!
         *   sets additional data for local elem FE.
         */
        virtual void
        set_local_fe_data(MAST::LocalElemFE& fe,
                          const libMesh::Elem& e) const;

        
        /*!
         *   clears the data structure of any stored values so that it can be 
         *   used for another element.
         */
        void clear();
        
        /*!
         *   calculates stress data for the element \p elem with the
         *   element solution \p sol.
         *   \p this object should have been cleared before calling this method.
         */
        void calculate_for_element(const libMesh::Elem& elem,
                                   const RealVectorX& sol);


        /*!
         *   calculates sensitivity of stress data for the element \p elem
         *   with the element solution \p sol.
         *   \p this object should have been cleared before calling this method.
         */
        void calculate_sensitivity_for_element(const libMesh::Elem& elem,
                                               const RealVectorX& sol,
                                               const RealVectorX& dsol,
                                               const MAST::FunctionBase& p);

        
        /*!
         *   @returns the thermal load for this element, if present in the
         *   the volume loads. A \p nullptr pointer is returned if no load is 
         *   present.
         */
        MAST::BoundaryConditionBase*
        get_thermal_load_for_elem(const libMesh::Elem& elem);

        
        /*!
         *   @returns the number of points for which stress-strain data is 
         *   stored for the given element.
         */
        unsigned int
        n_stress_strain_data_for_elem(const libMesh::Elem* e) const;

        
        
        /*!
         *   add the stress tensor associated with the qp. @returns a reference
         *   to \p Data.
         */
        MAST::StressStrainOutputBase::Data&
        add_stress_strain_at_qp_location(const libMesh::Elem*,
                                         const libMesh::Point& quadrature_pt,
                                         const libMesh::Point& physical_pt,
                                         const RealVectorX& strain,
                                         const RealVectorX& stress,
                                         Real JxW);
        
        
        /*!
         *    @returns the map of stress/strain data for all elems
         */
        const std::map<const libMesh::Elem*,
        std::vector<MAST::StressStrainOutputBase::Data*> >&
        get_stress_strain_data() const;

        
        /*!
         *    @returns the vector of stress/strain data for specified elem.
         */
        const std::vector<MAST::StressStrainOutputBase::Data*>&
        get_stress_strain_data_for_elem(const libMesh::Elem* e) const;
        
        
        
        /*!
         *   calculates and returns the von Mises p-norm functional for 
         *   all the elements that this object currently stores data for
         */
        Real
        von_Mises_p_norm_functional_for_all_elems(const Real p) const;

        
        /*!
         *   calculates and returns the sensitivity of von Mises p-norm
         *   functional for all the elements that this object currently 
         *   stores data for
         */
        Real
        von_Mises_p_norm_functional_sensitivity_for_all_elems
        (const Real p,
         const MAST::FunctionBase* f) const;

        
        /*!
         *   calculates and returns the derivative of von Mises p-norm
         *   functional wrt state vector for all the elements that
         *   this object currently stores data for
         */
        RealVectorX
        von_Mises_p_norm_functional_state_derivartive_for_all_elems(const Real p) const;

        
        
    protected:

        /*!
         *   \f$ p-\f$norm to be used for calculation of output stress function
         */
        Real _p_norm;
        
        
        /*!
         *    vector of stress with the associated location details
         */
        std::map<const libMesh::Elem*, std::vector<MAST::StressStrainOutputBase::Data*> >
        _stress_data;
        
    };
}

#endif // __mast__stress_output_base__
