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

#ifndef __mast__bending_operator__
#define __mast__bending_operator__

// C++ includes
#include <memory>

// MAST includes
#include "base/mast_data_types.h"

// libMesh includes
#include "libmesh/elem.h"
#include "libmesh/quadrature.h"
#include "libmesh/fe_base.h"



namespace MAST {
    
    // Forward declerations
    class StructuralElementBase;
    class FunctionBase;
    class SensitivityParameters;
    class FEMOperatorMatrix;
    class FEBase;
    template <typename ValType> class FieldFunction;
    
    
    enum BendingOperatorType {
        BERNOULLI,        // beam
        TIMOSHENKO,       // beam
        DKT,              // plate
        MINDLIN,          // plate
        DEFAULT_BENDING,
        NO_BENDING
    };
    
    
    class BendingOperator {
    public:
        
        BendingOperator(MAST::StructuralElementBase& elem);
        
        virtual ~BendingOperator();
        
        /*!
         *   returns true if this bending operator supports a transverse shear component
         */
        virtual bool include_transverse_shear_energy() const = 0;
        
        /*!
         *   calculate the transverse shear component for the element
         */
        virtual void
        calculate_transverse_shear_residual(bool request_jacobian,
                                            RealVectorX& local_f,
                                            RealMatrixX& local_jac)
        { libmesh_error(); }

        
        /*!
         *   calculate the transverse shear component for the element
         */
        virtual void
        calculate_transverse_shear_residual_sensitivity(const MAST::FunctionBase& p,
                                                        bool request_jacobian,
                                                        RealVectorX& local_f,
                                                        RealMatrixX& local_jac)
        { libmesh_error(); }


        /*!
         *   calculate the transverse shear component for the element
         */
        virtual void
        calculate_transverse_shear_residual_boundary_velocity
        (const MAST::FunctionBase& p,
         const unsigned int s,
         const MAST::FieldFunction<RealVectorX>& vel_f,
         bool request_jacobian,
         RealVectorX& local_f,
         RealMatrixX& local_jac)
        { libmesh_error(); }

        
    protected:
        
        /*!
         *   structural element associated with this
         */
        MAST::StructuralElementBase& _structural_elem;
        
        /*!
         *    element for which bending operator is created
         */
        const libMesh::Elem& _elem;
    };
    
    
    
    /*!
     *   Bending strain operator for 1D element
     */
    class BendingOperator1D:
    public MAST::BendingOperator {
        
    public:
        BendingOperator1D(MAST::StructuralElementBase& elem):
        MAST::BendingOperator(elem)
        { }
        
        /*!
         *   initialze the bending strain operator for the specified quadrature point
         */
        virtual void
        initialize_bending_strain_operator (const MAST::FEBase& fe,
                                            const unsigned int qp,
                                            MAST::FEMOperatorMatrix& Bmat_v,
                                            MAST::FEMOperatorMatrix& Bmat_w) = 0;

        /*!
         *   initialze the bending strain operator for the specified quadrature point
         */
        virtual void
        initialize_bending_strain_operator_for_yz (const MAST::FEBase& fe,
                                                   const unsigned int qp,
                                                   const Real y,
                                                   const Real z,
                                                   MAST::FEMOperatorMatrix& Bmat_v,
                                                   MAST::FEMOperatorMatrix& Bmat_w) = 0;
        
    };
    
    
    /*!
     *   Bending strain operator for 1D element
     */
    class BendingOperator2D: public MAST::BendingOperator {
    public:
        BendingOperator2D(MAST::StructuralElementBase& elem):
        MAST::BendingOperator(elem)
        { }
        
        /*!
         *   initialze the bending strain operator for the specified quadrature point
         */
        virtual void
        initialize_bending_strain_operator (const MAST::FEBase& fe,
                                            const unsigned int qp,
                                            MAST::FEMOperatorMatrix& Bmat) = 0;

        /*!
         *   initialze the bending strain operator for the specified quadrature point
         */
        virtual void
        initialize_bending_strain_operator_for_z (const MAST::FEBase& fe,
                                                  const unsigned int qp,
                                                  const Real z,
                                                  MAST::FEMOperatorMatrix& Bmat) = 0;
        
    };
    
    
    /*!
     *   builds a bending operator and returns it in a smart-pointer
     */
    std::unique_ptr<MAST::BendingOperator1D>
    build_bending_operator_1D(MAST::BendingOperatorType type,
                              MAST::StructuralElementBase& elem,
                              const std::vector<libMesh::Point>& pts);

    std::unique_ptr<MAST::BendingOperator2D>
    build_bending_operator_2D(MAST::BendingOperatorType type,
                              MAST::StructuralElementBase& elem,
                              const std::vector<libMesh::Point>& pts);

}

#endif

