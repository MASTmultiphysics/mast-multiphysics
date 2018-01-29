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

#ifndef __mast__output_function_base__
#define __mast__output_function_base__

// C++ includes
#include <set>
#include <memory>

// MAST includes
#include "base/mast_data_types.h"
#include "base/assembly_elem_operation.h"

// libMesh includes
#include "libmesh/elem.h"



namespace MAST {
    
    // Forward declerations
    class FunctionBase;
    class ElementBase;
    
    /*!
     *   This provides the base class for definitin of element level
     *   contribution of output quantity in an analysis.
     *   Apart from providing evaluation of the quantity
     *   \f$ q(X, p) \f$, the implementations are expected to evaluate
     *   the sensitivity \f$ \frac{\partial q}{\partial p} \f$ and the
     *   derivative \f$ \frac{\partial q}{\partial X} \f$. The latter quantity
     *   is used in an adjoint analysis.
     */
    class OutputAssemblyElemOperations:
    public MAST::AssemblyElemOperations {
        
    public:
        
        OutputAssemblyElemOperations();

        /*!
         *    virtual destructor
         */
        virtual ~OutputAssemblyElemOperations();
       
        
        /*!
         *   zeroes the output quantity values stored inside this object
         *   so that assembly process can begin.
         */
        virtual void zero() = 0;
        
        /*!
         *   @returns the output quantity value
         */
        virtual Real output() = 0;
        
        /*!
         *   @returns the output quantity sensitivity for parameter
         */
        virtual Real output_sensitivity(const MAST::FunctionBase& p) = 0;


        /*!
         *   @returns the output quantity derivative with respect to
         *   state vector
         */
        virtual RealVectorX output_derivative() = 0;

        
        /*!
         *    this is the abstract interface to be implemented by derived
         *    classes. This method calculates the quantity \f$ q(X, p) \f$.
         */
        virtual void evaluate() = 0;


        /*!
         *    this is the abstract interface to be implemented by derived
         *    classes. This method calculates the partial derivative of quantity
         *    \f[ \frac{\partial q(X, p)}{\partial p} \f]  with
         *    respect to parameter \f$ p \f$.
         */
        virtual void evaluate_sensitivity(const MAST::FunctionBase& p) = 0;

        
        /*!
         *    this is the abstract interface to be implemented by derived
         *    classes. This method calculates the quantity
         *    \f[ \frac{\partial q(X, p)}{\partial X} \f] for this
         *    output function.
         */
        virtual void evaluate_derivative() = 0;


        /*!
         *   The output function can be a boundary integrated quantity, volume
         *   integrated quantity or a combination of these two. The user
         *   must identify the subdomains on which the quantity will be defined
         *   using this method, or the \p set_participating_elements() method.
         */
        void set_participating_subdomains(const std::set<libMesh::subdomain_id_type>& sids);

        /*!
         *   sets the elements for which this object will evaluate and store
         *   the output data. This allows the user to specify a
         *   smaller subset of elements that will be grouped together in the
         *   output for constraint evaluation. If this method is
         *   not called, then the object will store data for all elements in
         *   the subdomain.
         */
        void set_participating_elements(const std::set<const libMesh::Elem*>& elems);

        
        /*!
         *   This will allow volume contribution from all elements.
         */
        void set_participating_elements_to_all();

        /*!
         *   The assembly will integration over boudnaries with ids specified in
         *   \p bids.
         */
        void set_participating_boundaries(const std::set<libMesh::boundary_id_type>& bids);

        /*!
         *    @returns the set of elements for which data will be stored. This
         *    is set using the \par set_participating_elements() method.
         */
        const std::set<const libMesh::Elem*>&
        get_participating_elements() const;
        

        /*!
         *    @returns the set of subdomain ids for which data will be stored.
         *    This is set using the \par set_participating_subdomains() method.
         */
        const std::set<libMesh::subdomain_id_type>&
        get_participating_subdomains();
        
        
        /*!
         *    @returns the set of boundary ids for which data will be stored.
         *    This is set using the \par set_participating_boundaries() method.
         */
        const std::set<libMesh::boundary_id_type>&
        get_participating_boundaries();

        
        /*!
         *    checks to see if the object has been told about the subset of
         *    elements and if the specified element is in the subset.
         */
        bool if_evaluate_for_element(const libMesh::Elem& elem) const;

        
    protected:
        
        /*!
         *   if true, evaluates on all elements.
         */
        bool _if_evaluate_on_all_elems;
        
        /*!
         *    set of elements for which the data will be stored. If this is
         *    empty, then data for all elements will be stored.
         */
        std::set<const libMesh::Elem*> _elem_subset;

        
        /*!
         *    set of subdomain ids for which data will be processed.
         */
        std::set<libMesh::subdomain_id_type> _sub_domain_ids;

        
        /*!
         *    set of bids for which data will be processed
         */
        std::set<libMesh::boundary_id_type> _bids;

    };
}

#endif // __mast__output_function_base__
