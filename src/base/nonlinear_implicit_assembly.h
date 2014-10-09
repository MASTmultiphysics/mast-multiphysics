
#ifndef __mast__nonlinear_implicit_assembly__
#define __mast__nonlinear_implicit_assembly__

// MAST includes
#include "base/assembly_base.h"

// libMesh includes
#include "libmesh/nonlinear_implicit_system.h"


namespace MAST {
    
    class NonlinearImplicitAssembly:
    public MAST::AssemblyBase,
    public libMesh::NonlinearImplicitSystem::ComputeResidualandJacobian,
    public libMesh::System::SensitivityAssembly {
    public:
        
        /*!
         *   constructor associates this assembly object with the system
         */
        NonlinearImplicitAssembly();
        
        
        /*!
         *   destructor resets the association of this assembly object with
         *   the system
         */
        virtual ~NonlinearImplicitAssembly();
        
        
        /*!
         *   attaches a system to this discipline, and vice-a-versa
         */
        virtual void
        attach_discipline_and_system(MAST::PhysicsDisciplineBase& discipline,
                                     MAST::SystemInitialization& system);
        
        
        /*!
         *   clears association with a system to this discipline, and vice-a-versa
         */
        virtual void
        clear_discipline_and_system( );

        
        /*!
         *    function that assembles the matrices and vectors quantities for
         *    nonlinear solution
         */
        virtual void
        residual_and_jacobian (const libMesh::NumericVector<Real>& X,
                               libMesh::NumericVector<Real>* R,
                               libMesh::SparseMatrix<Real>*  J,
                               libMesh::NonlinearImplicitSystem& S);
        
        /**
         * Assembly function.  This function will be called
         * to assemble the sensitivity of system residual prior to a solve and must
         * be provided by the user in a derived class. The method provides dR/dp_i
         * for \par i ^th parameter in the vector \par parameters.
         *
         * If the routine is not able to provide sensitivity for this parameter,
         * then it should return false, and the system will attempt to use
         * finite differencing.
         */
        virtual bool
        sensitivity_assemble (const libMesh::ParameterVector& parameters,
                              const unsigned int i,
                              libMesh::NumericVector<Real>& sensitivity_rhs);
        
    protected:
        
        /*!
         *   performs the element calculations over \par elem, and returns
         *   the element vector and matrix quantities in \par mat and
         *   \par vec, respectively. \par if_jac tells the method to also
         *   assemble the Jacobian, in addition to the residual vector.
         */
        virtual void _elem_calculations(MAST::ElementBase& elem,
                                        bool if_jac,
                                        RealVectorX& vec,
                                        RealMatrixX& mat) = 0;
        
        /*!
         *   performs the element sensitivity calculations over \par elem,
         *   and returns the element residual sensitivity in \par vec .
         */
        virtual void _elem_sensitivity_calculations(MAST::ElementBase& elem,
                                                    RealVectorX& vec) = 0;
    };
}


#endif //__mast__nonlinear_implicit_assembly__
