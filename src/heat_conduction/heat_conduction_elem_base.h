

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
                                  const MAST::ElementPropertyCardBase& p);
        
        
        virtual ~HeatConductionElementBase();
        
        
        /*!
         *   returns a constant reference to the finite element object
         */
        const MAST::ElementPropertyCardBase&
        elem_property()  {
            return _property;
        }
        
        
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
        
    protected:
        
                
        
        /*!
         *    Calculates the force vector and Jacobian due to surface pressure.
         *    this should be implemented for each element type
         */
        virtual bool surface_flux_residual(bool request_jacobian,
                                           RealVectorX& f,
                                           RealMatrixX& jac,
                                           MAST::BoundaryConditionBase& p,
                                           const libMesh::FEBase& fe);
        
        
        /*!
         *    Calculates the force vector and Jacobian due to surface pressure.
         */
        virtual bool
        surface_flux_residual_sensitivity(bool request_jacobian,
                                          RealVectorX& f,
                                          RealMatrixX& jac,
                                          MAST::BoundaryConditionBase& p,
                                          const libMesh::FEBase& fe);
        
        
        /*!
         *    Calculates the force vector and Jacobian due to surface pressure.
         *    this should be implemented for each element type
         */
        virtual bool surface_convection_residual(bool request_jacobian,
                                                 RealVectorX& f,
                                                 RealMatrixX& jac,
                                                 MAST::BoundaryConditionBase& p,
                                                 const libMesh::FEBase& fe);
        
        
        /*!
         *    Calculates the force vector and Jacobian due to surface pressure.
         */
        virtual bool
        surface_convection_residual_sensitivity(bool request_jacobian,
                                                RealVectorX& f,
                                                RealMatrixX& jac,
                                                MAST::BoundaryConditionBase& p,
                                                const libMesh::FEBase& fe);
        
        
        /*!
         *    Calculates the force vector and Jacobian due to surface pressure.
         *    this should be implemented for each element type
         */
        virtual bool surface_radiation_residual(bool request_jacobian,
                                                RealVectorX& f,
                                                RealMatrixX& jac,
                                                MAST::BoundaryConditionBase& p,
                                                const libMesh::FEBase& fe);
        
        
        /*!
         *    Calculates the force vector and Jacobian due to surface pressure.
         */
        virtual bool
        surface_radiation_residual_sensitivity(bool request_jacobian,
                                               RealVectorX& f,
                                               RealMatrixX& jac,
                                               MAST::BoundaryConditionBase& p,
                                               const libMesh::FEBase& fe);
        
        
        /*!
         *    Calculates the force vector and Jacobian due to surface pressure
         *    applied on the entire element domain. This is applicable for
         *    only 1D and 2D elements.
         */
        virtual bool volume_heat_source_residual(bool request_jacobian,
                                                 RealVectorX& f,
                                                 RealMatrixX& jac,
                                                 MAST::BoundaryConditionBase& p,
                                                 const libMesh::FEBase& fe);

        
        /*!
         *    Calculates the force vector and Jacobian due to surface pressure.
         *    this should be implemented for each element type
         */
        virtual bool
        volume_heat_source_residual_sensitivity(bool request_jacobian,
                                                RealVectorX& f,
                                                RealMatrixX& jac,
                                                MAST::BoundaryConditionBase& p,
                                                const libMesh::FEBase& fe);
        
        /*!
         *    When \p mass = false, initializes the FEM operator matrix to the
         *    shape functions as
         *    \f[  B = \left[ \begin{array}{c}
         *    {\bf N} \\ {\bf N} \\ {\bf N}
         *    \end{array} \right] \f]
         */
        void _initialize_mass_fem_operator(const unsigned int qp,
                                           MAST::FEMOperatorMatrix& Bmat);

        
        /*!
         *    For \p mass = true, the FEM operator matrix is initilized to
         *    the weak form of the Laplacian
         *    \f[  dB[0] = \frac{\partial {\bf N}}{\partial x} \f]
         *
         *    \f[  dB[1] = \frac{\partial {\bf N}}{\partial y} \f]
         *
         *    \f[  dB[2] = \frac{\partial {\bf N}}{\partial z} \f]
         */
        void _initialize_flux_fem_operator(const unsigned int qp,
                                           std::vector<MAST::FEMOperatorMatrix>& dBmat);

        /*!
         *   element property
         */
        const MAST::ElementPropertyCardBase& _property;
        
    };
    
}


#endif // __mast__heat_conduction_elem_base__
