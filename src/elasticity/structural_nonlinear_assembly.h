
#ifndef __mast__structural_nonlinear_assembly__
#define __mast__structural_nonlinear_assembly__

// MAST includes
#include "base/nonlinear_implicit_assembly.h"


namespace MAST {
    
    
    class StructuralNonlinearAssembly:
    public MAST::NonlinearImplicitAssembly {
        
    public:
        
        /*!
         *   constructor associates this assembly object with the system
         */
        StructuralNonlinearAssembly();
        
        
        /*!
         *   destructor resets the association of this assembly object with
         *   the system
         */
        virtual ~StructuralNonlinearAssembly();
        

        /*!
         *    function that assembles the matrices and vectors quantities for
         *    nonlinear solution
         */
        virtual void
        residual_and_jacobian (const libMesh::NumericVector<Real>& X,
                               libMesh::NumericVector<Real>* R,
                               libMesh::SparseMatrix<Real>*  J,
                               libMesh::NonlinearImplicitSystem& S);

        /*!
         *   asks the system to update the nonlinear incompatible mode solution
         */
        void update_incompatible_solution(libMesh::NumericVector<Real>& X,
                                          libMesh::NumericVector<Real>& dX);
        
    protected:
        
        /*!
         *   @returns a smart-pointer to a newly created element for
         *   calculation of element quantities.
         */
        virtual std::auto_ptr<MAST::ElementBase>
        _build_elem(const libMesh::Elem& elem);
        
        /*!
         *   performs the element calculations over \par elem, and returns
         *   the element vector and matrix quantities in \par mat and
         *   \par vec, respectively. \par if_jac tells the method to also
         *   assemble the Jacobian, in addition to the residual vector.
         */
        virtual void _elem_calculations(MAST::ElementBase& elem,
                                        bool if_jac,
                                        RealVectorX& vec,
                                        RealMatrixX& mat);
        
        /*!
         *   performs the element sensitivity calculations over \par elem,
         *   and returns the element residual sensitivity in \par vec .
         */
        virtual void _elem_sensitivity_calculations(MAST::ElementBase& elem,
                                                    RealVectorX& vec);
        
        /*!
         *   map of local incompatible mode solution per 3D elements
         */
        std::map<const libMesh::Elem*, RealVectorX> _incompatible_sol;
    };
}


#endif // __mast__structural_nonlinear_assembly__
