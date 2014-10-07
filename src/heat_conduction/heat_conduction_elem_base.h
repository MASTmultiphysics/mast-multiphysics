

#ifndef __mast__heat_conduction_elem_base__
#define __mast__heat_conduction_elem_base__

// MAST includes
#include "base/elem_base.h"

namespace MAST {
    
    // Forward declerations
    class ElementPropertyCardBase;
    class LocalElemBase;
    class BoundaryConditionBase;
    class FEMOperatorMatrix;
    
    class HeatConductionElementBase:
    public MAST::ElementBase
    {
    public:
        /*!
         *   Constructor
         */
        HeatConductionElementBase(MAST::SystemInitialization& sys,
                                  const libMesh::Elem& elem,
                                  MAST::ElementPropertyCardBase& p);
        
        
        virtual ~HeatConductionElementBase();
        
        
        /*!
         *   returns a constant reference to the finite element object
         */
        const MAST::ElementPropertyCardBase&
        elem_property()  {
            return _property;
        }
        
        
        /*!
         *   returns a constant reference to the element in local coordinate system
         */
        virtual const libMesh::Elem&
        get_elem_for_quadrature() const = 0;
        
        
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
                                       RealMatrixX& jac,
                                       bool if_ignore_ho_jac) = 0;
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
        
    protected:
        
                
        
        /*!
         *    Calculates the force vector and Jacobian due to surface pressure.
         *    this should be implemented for each element type
         */
        virtual bool surface_flux_residual(bool request_jacobian,
                                           RealVectorX& f,
                                           RealMatrixX& jac,
                                           const unsigned int side,
                                           MAST::BoundaryConditionBase& p);
        
        
        /*!
         *    Calculates the force vector and Jacobian due to surface pressure
         *    applied on the entire element domain. This is applicable for
         *    only 1D and 2D elements.
         */
        virtual bool surface_flux_residual(bool request_jacobian,
                                           RealVectorX& f,
                                           RealMatrixX& jac,
                                           MAST::BoundaryConditionBase& p);
        
        
        
        /*!
         *    Calculates the force vector and Jacobian due to surface pressure.
         */
        virtual bool surface_flux_residual_sensitivity(bool request_jacobian,
                                                       RealVectorX& f,
                                                       RealMatrixX& jac,
                                                       const unsigned int side,
                                                       MAST::BoundaryConditionBase& p);
        
        
        /*!
         *    Calculates the force vector and Jacobian due to surface pressure.
         *    this should be implemented for each element type
         */
        virtual bool surface_flux_residual_sensitivity(bool request_jacobian,
                                                       RealVectorX& f,
                                                       RealMatrixX& jac,
                                                       MAST::BoundaryConditionBase& p);

        /*!
         *    Calculates the force vector and Jacobian due to surface pressure.
         *    this should be implemented for each element type
         */
        virtual bool surface_convection_residual(bool request_jacobian,
                                                 RealVectorX& f,
                                                 RealMatrixX& jac,
                                                 const unsigned int side,
                                                 MAST::BoundaryConditionBase& p);
        
        
        /*!
         *    Calculates the force vector and Jacobian due to surface pressure
         *    applied on the entire element domain. This is applicable for
         *    only 1D and 2D elements.
         */
        virtual bool surface_convection_residual(bool request_jacobian,
                                                 RealVectorX& f,
                                                 RealMatrixX& jac,
                                                 MAST::BoundaryConditionBase& p);
        
        
        
        /*!
         *    Calculates the force vector and Jacobian due to surface pressure.
         */
        virtual bool surface_convection_residual_sensitivity(bool request_jacobian,
                                                             RealVectorX& f,
                                                             RealMatrixX& jac,
                                                             const unsigned int side,
                                                             MAST::BoundaryConditionBase& p);
        
        
        /*!
         *    Calculates the force vector and Jacobian due to surface pressure.
         *    this should be implemented for each element type
         */
        virtual bool surface_convection_residual_sensitivity(bool request_jacobian,
                                                             RealVectorX& f,
                                                             RealMatrixX& jac,
                                                             MAST::BoundaryConditionBase& p);

        
        /*!
         *    Calculates the force vector and Jacobian due to surface pressure.
         *    this should be implemented for each element type
         */
        virtual bool surface_radiation_residual(bool request_jacobian,
                                                RealVectorX& f,
                                                RealMatrixX& jac,
                                                const unsigned int side,
                                                MAST::BoundaryConditionBase& p);
        
        
        /*!
         *    Calculates the force vector and Jacobian due to surface pressure
         *    applied on the entire element domain. This is applicable for
         *    only 1D and 2D elements.
         */
        virtual bool surface_radiation_residual(bool request_jacobian,
                                                RealVectorX& f,
                                                RealMatrixX& jac,
                                                MAST::BoundaryConditionBase& p);
        
        
        
        /*!
         *    Calculates the force vector and Jacobian due to surface pressure.
         */
        virtual bool surface_radiation_residual_sensitivity(bool request_jacobian,
                                                            RealVectorX& f,
                                                            RealMatrixX& jac,
                                                            const unsigned int side,
                                                            MAST::BoundaryConditionBase& p);
        
        
        /*!
         *    Calculates the force vector and Jacobian due to surface pressure.
         *    this should be implemented for each element type
         */
        virtual bool surface_radiation_residual_sensitivity(bool request_jacobian,
                                                            RealVectorX& f,
                                                            RealMatrixX& jac,
                                                            MAST::BoundaryConditionBase& p);
        

        /*!
         *    Calculates the force vector and Jacobian due to surface pressure
         *    applied on the entire element domain. This is applicable for
         *    only 1D and 2D elements.
         */
        virtual bool volume_heat_source_residual(bool request_jacobian,
                                                 RealVectorX& f,
                                                 RealMatrixX& jac,
                                                 MAST::BoundaryConditionBase& p);

        
        /*!
         *    Calculates the force vector and Jacobian due to surface pressure.
         *    this should be implemented for each element type
         */
        virtual bool volume_heat_source_residual_sensitivity(bool request_jacobian,
                                                             RealVectorX& f,
                                                             RealMatrixX& jac,
                                                             MAST::BoundaryConditionBase& p);

        /*!
         *    When \p mass = false, initializes the FEM operator matrix to the
         *    shape functions as
         *    \f[  B = \left[ \begin{array}{c}
         *    {\bf N} \\ {\bf N} \\ {\bf N}
         *    \end{array} \right] \f]
         *
         *    For \p mass = true, the FEM operator matrix is initilized to
         *    the weak form of the Laplacian
         *    \f[  B = \left[ \begin{array}{c}
         *    \frac{\partial {\bf N}}{\partial x} \\
         *    \frac{\partial {\bf N}}{\partial x} \\
         *    \frac{\partial {\bf N}}{\partial x}
         *    \end{array} \right] \f]
         */
        void _initialize_fem_operator(const unsigned int qp,
                                      const bool mass,
                                      MAST::FEMOperatorMatrix& Bmat);
        
        /*!
         *   element property
         */
        MAST::ElementPropertyCardBase& _property;
        
    };
    
}


#endif // __mast__heat_conduction_elem_base__
