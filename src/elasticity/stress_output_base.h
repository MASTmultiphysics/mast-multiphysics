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

#ifndef __mast__stress_output_base__
#define __mast__stress_output_base__

// C++ includes
#include <map>

// MAST includes
#include "base/output_function_base.h"
#include "base/mast_data_types.h"


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
     *    stresses in the x, y, z directions, sigma_xx, sigma_yy, sigma_zz, and
     *    the next three components are the shear stresses
     *    tau_xy, tau_yz, tau_xz.
     *
     *    The first three components of strain are the direct
     *    strains in the x, y, z directions, epsilon_xx, epsilon_yy, epsilon_zz,
     *    and the next three components are the engineering shear strains
     *    gamma_xy, gamma_yz, gamma_xz.

     */
    class StressStrainOutputBase:
    public MAST::OutputFunctionBase {
        
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
                 const libMesh::Point& xyz);
 
            
            
            
            /*!
             *   @returns stress
             */
            const RealVectorX stress() const;

            
            /*!
             *   @returns strain
             */
            const RealVectorX strain() const;
            
            
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

            
            libMesh::Point            _qp;
            
            
            libMesh::Point            _xyz;
        };
        

        
        /*!
         *    default constructor
         */
        StressStrainOutputBase();
        
        virtual ~StressStrainOutputBase();
        

        /*!
         *   clears the data structure of any stored values so that it can be 
         *   used for another element.
         */
        void clear();
        
        
        /*!
         *   add the stress tensor associated with the qp. @returns a reference
         *   to \p Data.
         */
        MAST::StressStrainOutputBase::Data&
        add_stress_strain_at_qp_location(const libMesh::Point& quadrature_pt,
                                         const libMesh::Point& physical_pt,
                                         const RealVectorX& strain,
                                         const RealVectorX& stress);
        
        
        /*!
         *    @returns the vector of stress/strain data
         */
        const std::vector<MAST::StressStrainOutputBase::Data*>&
        get_stress_strain_data() const;
        
    protected:

        
        /*!
         *    vector of stress with the associated location details
         */
        std::vector<MAST::StressStrainOutputBase::Data*> _stress_data;

    };
}

#endif // __mast__stress_output_base__
