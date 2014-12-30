/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2014  Manav Bhatia
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

#ifndef __mast__solid_1d_section_element_property_card__
#define __mast__solid_1d_section_element_property_card__


// MAST includes
#include "property_cards/element_property_card_1D.h"

namespace MAST {
    
    
    class Solid1DSectionElementPropertyCard :
    public MAST::ElementPropertyCard1D {
        
    public:
        
        Solid1DSectionElementPropertyCard():
        MAST::ElementPropertyCard1D(),
        _initialized(false),
        _material(NULL)
        { }
        
        
        /*!
         *   virtual destructor
         */
        virtual ~Solid1DSectionElementPropertyCard() { }
        
        
        virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
        stiffness_A_matrix(const MAST::ElementBase& e);
        
        virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
        stiffness_B_matrix(const MAST::ElementBase& e);
        
        virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
        stiffness_D_matrix(const MAST::ElementBase& e);
        
        virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
        damping_matrix(const MAST::ElementBase& e);
        
        virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
        inertia_matrix(const MAST::ElementBase& e);
        
        virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
        thermal_expansion_A_matrix(const MAST::ElementBase& e);
        
        virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
        thermal_expansion_B_matrix(const MAST::ElementBase& e);
        
        virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
        transverse_shear_stiffness_matrix(const MAST::ElementBase& e);
        
        virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
        prestress_A_matrix(const MAST::ElementBase& e);
        
        virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
        prestress_B_matrix(const MAST::ElementBase& e);

        
        virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
        thermal_conductance_matrix(const MAST::ElementBase& e);
        

        virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
        thermal_capacitance_matrix(const MAST::ElementBase& e);

        
        /*!
         *    sets the material card
         */
        void set_material(MAST::MaterialPropertyCardBase& mat) {
            _material = &mat;
        }
        
        
        /*!
         *   return true if the property is isotropic
         */
        virtual bool if_isotropic() const {
            return true;
        }
        
        
        /*!
         *    returns a reference to the material
         */
        virtual const MAST::MaterialPropertyCardBase& get_material() const {
            libmesh_assert(_material); // make sure it has already been set
            return *_material;
        }
        
        /*!
         *   @returns reference to the section area function
         */
        MAST::FieldFunction<Real>& A();

        /*!
         *   @returns reference to the function that calculates the 
         *   section torsional constant
         */
        MAST::FieldFunction<Real>& J();

        /*!
         *   @returns reference to the section polar moment function
         */
        MAST::FieldFunction<Real>& Ip();

        /*!
         *   @returns reference to the function that calculates
         *   section area moment about y-axis
         */
        MAST::FieldFunction<Real>& Ay();

        /*!
         *   @returns reference to the function that calculates
         *   section area moment about z-axis
         */
        MAST::FieldFunction<Real>& Az();

        /*!
         *   @returns reference to the function that calculates the 
         *   section area moment of inertia
         */
        MAST::FieldFunction<RealMatrixX>& I();

        /*!
         *  returns true if the property card depends on the function \p f
         */
        virtual bool depends_on(const MAST::FunctionBase& f) const;

        
        virtual void init();
        
    protected:

        bool _initialized;
        
        /*!
         *   material property card
         */
        MAST::MaterialPropertyCardBase *_material;
        
        std::auto_ptr<MAST::FieldFunction<Real> > _A;
        
        std::auto_ptr<MAST::FieldFunction<Real> > _J;

        std::auto_ptr<MAST::FieldFunction<Real> > _Ip;

        std::auto_ptr<MAST::FieldFunction<Real> > _Ay;
        
        std::auto_ptr<MAST::FieldFunction<Real> > _Az;
        
        std::auto_ptr<MAST::FieldFunction<RealMatrixX> > _AI;
        
        std::auto_ptr<MAST::FieldFunction<RealMatrixX> > _stiff_A;

        std::auto_ptr<MAST::FieldFunction<RealMatrixX> > _stiff_B;

        std::auto_ptr<MAST::FieldFunction<RealMatrixX> > _stiff_D;

        std::auto_ptr<MAST::FieldFunction<RealMatrixX> > _damp;

        std::auto_ptr<MAST::FieldFunction<RealMatrixX> > _inertia;

        std::auto_ptr<MAST::FieldFunction<RealMatrixX> > _thermal_A;

        std::auto_ptr<MAST::FieldFunction<RealMatrixX> > _thermal_B;

        std::auto_ptr<MAST::FieldFunction<RealMatrixX> > _transverse_shear;

    };
    
}



#endif // __mast__solid_1d_section_element_property_card__
