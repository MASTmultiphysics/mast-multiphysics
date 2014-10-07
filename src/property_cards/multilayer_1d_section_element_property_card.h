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

#ifndef __mast__multilayer_1d_section_element_property_card__
#define __mast__multilayer_1d_section_element_property_card__


// MAST includes
#include "property_cards/element_property_card_1D.h"


namespace MAST {

    // Forward declerations
    class Solid1DSectionElementPropertyCard;
    
    
    class Multilayer1DSectionElementPropertyCard :
    public MAST::ElementPropertyCard1D {
        
    public:
        
        Multilayer1DSectionElementPropertyCard(unsigned int pid):
        MAST::ElementPropertyCard1D(pid)
        { }
        
        
        /*!
         *   virtual destructor
         */
        virtual ~Multilayer1DSectionElementPropertyCard();
        
        
        /*!
         *    sets the layers of this section.
         *    \p base is used reference for calculation of offset.
         *    base = -1 implies section lower thickness,
         *    base implies section mid-point
         *    base = +1 implies section upper thickness.
         */
        void set_layers(const Real base,
                        std::vector<MAST::Solid1DSectionElementPropertyCard*>& layers);
        

        /*!
         *    returns the layers of this section
         */
        const std::vector<MAST::Solid1DSectionElementPropertyCard*>& get_layers() const;
        
        
        /*!
         *   return true if the property is isotropic
         */
        virtual bool if_isotropic() const;
        

        /*!
         *  returns true if the property card depends on the function \p f
         */
        virtual bool depends_on(const MAST::FunctionBase& f) const;
        
        
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

        
    protected:

        std::vector<MAST::FieldFunction<Real>*> _layer_offsets;
        
        /*!
         *   vector of thickness function for each layer
         */
        std::vector<MAST::Solid1DSectionElementPropertyCard*> _layers;
    };
    
}




#endif // __mast__multilayer_1d_section_element_property_card__
