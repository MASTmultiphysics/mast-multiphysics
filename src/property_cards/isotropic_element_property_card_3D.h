/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2020  Manav Bhatia and MAST authors
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

#ifndef __mast__isotropic_element_property_card_3D__
#define __mast__isotropic_element_property_card_3D__

// MAST includes
#include "property_cards/element_property_card_base.h"



namespace MAST
{
    
    class IsotropicElementPropertyCard3D:
    public MAST::ElementPropertyCardBase {
        
    public:
        IsotropicElementPropertyCard3D():
        MAST::ElementPropertyCardBase(),
        _material(nullptr)
        { }
        
        /*!
         *   virtual destructor
         */
        virtual ~IsotropicElementPropertyCard3D() { }
        
        
        /*!
         *   dimension of the element for which this property is defined
         */
        virtual unsigned int dim() const {
            return 3;
        }

        /*!
         *   return true if the property is isotropic
         */
        virtual bool if_isotropic() const {
            return true;
        }

        
        /*!
         *    sets the material card
         */
        virtual void set_material(MAST::MaterialPropertyCardBase& mat) {
            _material = &mat;
        }
        
        
        /*!
         *    returns a reference to the material
         */
        const MAST::MaterialPropertyCardBase& get_material() const {
            libmesh_assert(_material); // make sure it has already been set
            return *_material;
        }
        

        /*!
         *   returns the bending model to be used for the element. Should be
         *   reimplemented in the derived classes
         */
        virtual MAST::BendingOperatorType
        bending_model(const MAST::GeomElem& elem) const {
            libmesh_assert(false);
        }
        
        /*!
         *    returns the extra quadrature order (on top of the system) that
         *    this element should use. By default this is zero, and can be
         *    changed by the derived classes
         */
        virtual int extra_quadrature_order(const MAST::GeomElem& elem) const {
            libmesh_assert(false);
            return 0;
        }

        /*!
         *  returns true if the property card depends on the function \p f
         */
        virtual bool depends_on(const MAST::FunctionBase& f) const;

        
        virtual std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
        stiffness_A_matrix(const MAST::ElementBase& e) const;
        
        virtual std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
        stiffness_A_matrix() const;
        
        virtual std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
        stiffness_B_matrix(const MAST::ElementBase& e) const;
        
        virtual std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
        stiffness_B_matrix() const;
        
        virtual std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
        stiffness_D_matrix(const MAST::ElementBase& e) const;
        
        virtual std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
        stiffness_D_matrix() const;
        
        virtual std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
        damping_matrix(const MAST::ElementBase& e) const;
        
        virtual std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
        inertia_matrix(const MAST::ElementBase& e) const;
        
        virtual std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
        inertia_matrix() const;
        
        virtual std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
        thermal_expansion_A_matrix(const MAST::ElementBase& e) const;
        
        virtual std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
        thermal_expansion_A_matrix() const;
        
        virtual std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
        thermal_expansion_B_matrix(const MAST::ElementBase& e) const;
        
        virtual std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
        thermal_expansion_B_matrix() const;
        
        virtual std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
        transverse_shear_stiffness_matrix(const MAST::ElementBase& e) const;
        
        virtual std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
        prestress_A_matrix( MAST::ElementBase& e) const;
        
        virtual std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
        prestress_B_matrix( MAST::ElementBase& e) const;
        
        virtual std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
        thermal_conductance_matrix(const MAST::ElementBase& e) const;
        
        virtual std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
        thermal_conductance_matrix() const;
        
        virtual std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
        thermal_capacitance_matrix(const MAST::ElementBase& e) const;
        
        virtual std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
        thermal_capacitance_matrix() const;
        
        virtual const MAST::FieldFunction<Real>*
        section(const MAST::ElementBase& e) const {

            return nullptr;
        }

    protected:
        
        /*!
         *    pointer to the material property card
         */
        MAST::MaterialPropertyCardBase* _material;
    };
}


#endif // __mast__isotropic_element_property_card_3D__
