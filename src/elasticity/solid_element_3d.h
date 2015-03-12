//
//  solid_element_3d.h
//  mast
//
//  Created by Manav Bhatia on 2/17/15.
//  Copyright (c) 2015 MAST. All rights reserved.
//

#ifndef __mast__solid_element_3d__
#define __mast__solid_element_3d__

/// MAST includes
#include "elasticity/structural_element_base.h"

// Forward declerations
class FEMOperatorMatrix;

namespace MAST {
    
    // Forward declerations
    class BoundaryConditionBase;
    
    
    class StructuralElement3D:
    public MAST::StructuralElementBase {
        
    public:
        StructuralElement3D(MAST::SystemInitialization& sys,
                            const libMesh::Elem& elem,
                            const MAST::ElementPropertyCardBase& p):
        StructuralElementBase(sys, elem, p) {
            
        }
        
        /*!
         *    Calculates the internal residual vector and Jacobian due to
         *    strain energy
         */
        virtual bool internal_residual(bool request_jacobian,
                                       RealVectorX& f,
                                       RealMatrixX& jac,
                                       bool if_ignore_ho_jac);
        
        
        /*!
         *    Calculates the sensitivity of internal residual vector and
         *    Jacobian due to strain energy
         */
        virtual bool internal_residual_sensitivity(bool request_jacobian,
                                                   RealVectorX& f,
                                                   RealMatrixX& jac,
                                                   bool if_ignore_ho_jac);
        
        /*!
         *   calculates d[J]/d{x} . d{x}/dp
         */
        virtual bool
        internal_residual_jac_dot_state_sensitivity (RealMatrixX& jac) {
            libmesh_assert(false); // to be implemented for 3D elements
        }
        
        /*!
         *    Calculates the prestress residual vector and Jacobian
         */
        virtual bool prestress_residual (bool request_jacobian,
                                         RealVectorX& f,
                                         RealMatrixX& jac);
        
        
        /*!
         *    Calculates the sensitivity prestress residual vector and Jacobian
         */
        virtual bool prestress_residual_sensitivity (bool request_jacobian,
                                                     RealVectorX& f,
                                                     RealMatrixX& jac);
        
        
    protected:
        
        /*!
         *    Calculates the residual vector and Jacobian due to thermal stresses
         */
        virtual bool thermal_residual(bool request_jacobian,
                                      RealVectorX& f,
                                      RealMatrixX& jac,
                                      MAST::BoundaryConditionBase& p);
        
        /*!
         *    Calculates the sensitivity of residual vector and Jacobian due to
         *    thermal stresses
         */
        virtual bool thermal_residual_sensitivity(bool request_jacobian,
                                                  RealVectorX& f,
                                                  RealMatrixX& jac,
                                                  MAST::BoundaryConditionBase& p);
        
        /*!
         *   initialize strain operator matrix
         */
        void initialize_strain_operator (const unsigned int qp,
                                         FEMOperatorMatrix& Bmat);
    };
}

#endif // __mast__solid_element_3d__
