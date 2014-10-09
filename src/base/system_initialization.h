
#ifndef __mast__system_initialization__
#define __mast__system_initialization__


// libMesh includes
#include "libmesh/system.h"


namespace MAST {
    
    class SystemInitialization {
    public:
        /*!
         *   initialize the variables in the provided system \par sys
         *   of \par order and \par family. Uses \par prefix for
         *   all variables name.
         */
        SystemInitialization (libMesh::System& sys,
                              const std::string& prefix);
        
        /*!
         *   virtual destructor
         */
        virtual ~SystemInitialization();
        
        
        /*!
         *   @returns the number of variables in this system
         */
        unsigned int n_vars() const {
            return _system.n_vars();
        }
        
        /*!
         *   @returns the FEType object for variable \par i,
         *   that defines the finite element family and order.
         */
        const libMesh::FEType&
        fetype(unsigned int i) const {
            return _system.variable_type(i);
        }
        
        
        /*!
         *  @returns a reference to the system for which the variables
         *  are initialized.
         */
        libMesh::System& system() {
            return _system;
        }
        
        /*!
         *  @returns a constant reference to the system for which the variables
         *  are initialized.
         */
        const libMesh::System& system() const {
            return _system;
        }
        
        /*!
         *  @returns a constant reference to the vector of variable IDs.
         */
        const std::vector<unsigned int> vars() const {
            return _vars;
        }
        
        /*!
         *  @returns a constant reference to the prefix used for all
         *  variables.
         */
        const std::string& prefix() const {
            return _prefix;
        }
        
    protected:
        
        libMesh::System& _system;
        
        std::vector<unsigned int> _vars;
        
        std::string _prefix;
    };
}


#endif //__mast__system_initialization__
