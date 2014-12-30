
#ifndef __mast__mesh_field_function__
#define __mast__mesh_field_function__

// MAST includes
#include "base/field_function_base.h"


// libMesh includes
#include "libmesh/numeric_vector.h"
#include "libmesh/mesh_function.h"



namespace MAST {

    // Forward declerations
    class SystemInitialization;
    
    
    /*!
     *    This provides a wrapper FieldFunction compatible class that
     *    interpolates the solution using libMesh's MeshFunction class.
     */
    template <typename ValType>
    class MeshFieldFunction:
    public MAST::FieldFunction<ValType> {
        
    public:
        /*!
         *   constructor
         */
        MeshFieldFunction(const std::string& nm);
        

        /*!
         *   copy constructor
         */
        MeshFieldFunction(const MAST::MeshFieldFunction<ValType>& f);

        
        /*!
         *   destructor
         */
        ~MeshFieldFunction();

       
        /*!
         *   @returns a clone of the function
         */
        virtual std::auto_ptr<MAST::FieldFunction<ValType> > clone() const;
        
        /*!
         *    calculates the value of the function at the specified point,
         *    \par p, and time, \par t, and returns it in \p v.
         */
        virtual void operator() (const libMesh::Point& p,
                                 const Real t,
                                 ValType& v) const;
        
        
        /*!
         *    calculates the value of the function at the specified point,
         *    \par p, and time, \par t, and returns it in \p v.
         */
        virtual void derivative (const MAST::DerivativeType d,
                                 const MAST::FunctionBase& f,
                                 const libMesh::Point& p,
                                 const Real t,
                                 ValType& v) const;
        
        
        /*!
         *   initializes the data structures to perform the interpolation 
         *   function on the given \par system and \par sol.
         */
        void init_for_system_and_solution(MAST::SystemInitialization& sys,
                                          const libMesh::NumericVector<Real>& sol);

        
        /*!
         *    When a mesh field function is attached to an assembly routine
         *    during system assembly, then the current solution can be 
         *    provided by the element quadrature point update. This method
         *    allows the element to provide the solution and in turn
         *    override the mesh function evaluation.
         */
        void set_element_quadrature_point_solution(RealVectorX& sol);

        
        /*!
         *    clears the quadrature point solution provided by the 
         *    corresponding set method above.
         */
        void clear_element_quadrature_point_solution();

        
        /*!
         *   clear the solution of
         */
        void clear();

    protected:

        /*!
         *  flag is set to true when the quadrature point solution is 
         *  provided by an element
         */
        bool _use_qp_sol;
        
        
        /*!
         *   quadrature point solution of the element
         */
        RealVectorX _qp_sol;
        
        /*!
         *  current system for which solution is to be interpolated
         */
        MAST::SystemInitialization* _system;
        
        /*!
         *   current solution that is going to be interpolated
         */
        libMesh::NumericVector<Real>* _sol;
        
        /*!
         *   the MeshFunction object that performs the interpolation
         */
        libMesh::FunctionBase<Real>* _mesh_function;
    };
}

#endif // __mast__mesh_field_function__

