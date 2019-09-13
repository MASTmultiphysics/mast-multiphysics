/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2019  Manav Bhatia
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

/*
 * Below is a description of the cross sections that are implemented and the 
 * number of parameters which are used to define them.
 * 
 * Cross Section                # Parameters
 * --------------------         ------------
 * Square                       1
 * ROD                          1
 * Equilaterial Triangle        1
 * Regular Hexagon              1
 * BAR                          2
 * TUBE                         2
 * Ellipse                      2
 * Square Box                   2
 * BOX                          4   TODO
 * CHAN                         4   TODO
 * CROSS                        4   TODO
 * H                            4   TODO
 * L                            4   TODO
 * HAT (Open)                   4   TODO
 * I1                           4   TODO
 * T                            4   TODO
 * Z                            4   TODO
 * PI                           5   TODO
 * HAT1 (Closed)                5   TODO
 * BOX1                         6   TODO
 * I                            6   TODO
 * GBOX                         6   TODO
 * Arbitrary I-Beam             8   TODO
 * Arbitrary Section            Variable^
 * 
 *  ^   The arbitrary cross section is define via either specifying section 
 *      properties directly, via a matrix of vertices, or via a 2D 
 *      libMesh:MeshBase object which represents the cross section.
 * 
 * Note that thin-walled typtically means the thickness of the wall is less 
 * than 1/10 of the radius (1/20 of the diameter) of the circular tube. For 
 * rectangular tubes, define this as less than 1/20 of the shorest side.
 * 
 * References
 * ----------
 * https://structx.com/geometric_properties.html
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
        _material(nullptr)
        { }
        
        
        /*!
         *   virtual destructor
         */
        virtual ~Solid1DSectionElementPropertyCard() { }
        
        
        virtual std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
        stiffness_A_matrix(const MAST::ElementBase& e) const;
        
        virtual std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
        stiffness_B_matrix(const MAST::ElementBase& e) const;
        
        virtual std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
        stiffness_D_matrix(const MAST::ElementBase& e) const;
        
        virtual std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
        damping_matrix(const MAST::ElementBase& e) const;
        
        virtual std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
        inertia_matrix(const MAST::ElementBase& e) const;
        
        virtual std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
        thermal_expansion_A_matrix(const MAST::ElementBase& e) const;
        
        virtual std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
        thermal_expansion_B_matrix(const MAST::ElementBase& e) const;
        
        virtual std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
        transverse_shear_stiffness_matrix(const MAST::ElementBase& e) const;
        
        virtual std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
        prestress_A_matrix(MAST::ElementBase& e) const;
        
        virtual std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
        prestress_B_matrix(MAST::ElementBase& e) const;

        
        virtual std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
        thermal_conductance_matrix(const MAST::ElementBase& e) const;
        

        virtual std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
        thermal_capacitance_matrix(const MAST::ElementBase& e) const;

        virtual const MAST::FieldFunction<Real>&
        section(const MAST::ElementBase& e) const;

        /*!
         *    sets the material card
         */
        virtual void set_material(MAST::MaterialPropertyCardBase& mat) {
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
         *   @returns a constant reference to the section area function
         */
        virtual const MAST::FieldFunction<Real>& A() const;

        
        /*!
         *   @returns reference to the section area function
         */
        virtual MAST::FieldFunction<Real>& A();

        /*!
         *   @returns constant reference to the function that calculates the
         *   section torsional constant
         */
        virtual const MAST::FieldFunction<Real>& J() const;

        /*!
         *   @returns reference to the function that calculates the 
         *   section torsional constant
         */
        virtual MAST::FieldFunction<Real>& J();

        /*!
         *   @returns reference to the section polar moment function
         */
        virtual const MAST::FieldFunction<Real>& Ip() const;

        /*!
         *   @returns reference to the section polar moment function
         */
        virtual MAST::FieldFunction<Real>& Ip();

        /*!
         *   @returns constant reference to the function that calculates
         *   section area moment about y-axis
         */
        virtual const MAST::FieldFunction<Real>& Ay() const;

        /*!
         *   @returns reference to the function that calculates
         *   section area moment about y-axis
         */
        virtual MAST::FieldFunction<Real>& Ay();

        /*!
         *   @returns constant reference to the function that calculates
         *   section area moment about z-axis
         */
        virtual const MAST::FieldFunction<Real>& Az() const;

        /*!
         *   @returns reference to the function that calculates
         *   section area moment about z-axis
         */
        virtual MAST::FieldFunction<Real>& Az();

        /*!
         *   @returns constant reference to the function that calculates the
         *   section area moment of inertia
         */
        virtual const MAST::FieldFunction<RealMatrixX>& I() const;

        
        /*!
         *   @returns reference to the function that calculates the 
         *   section area moment of inertia
         */
        virtual MAST::FieldFunction<RealMatrixX>& I();

        /*!
         *  returns true if the property card depends on the function \p f
         */
        virtual bool depends_on(const MAST::FunctionBase& f) const;

        
        virtual void clear();
        
        virtual void init();
        
    protected:

        bool _initialized;
        
        /*!
         *   material property card
         */
        MAST::MaterialPropertyCardBase *_material;
        
        std::unique_ptr<MAST::FieldFunction<Real> > _A;
        
        std::unique_ptr<MAST::FieldFunction<Real> > _J;

        std::unique_ptr<MAST::FieldFunction<Real> > _Ip;

        std::unique_ptr<MAST::FieldFunction<Real> > _Ay;
        
        std::unique_ptr<MAST::FieldFunction<Real> > _Az;
        
        std::unique_ptr<MAST::FieldFunction<RealMatrixX> > _AI;
        
        std::unique_ptr<MAST::FieldFunction<RealMatrixX> > _stiff_A;

        std::unique_ptr<MAST::FieldFunction<RealMatrixX> > _stiff_B;

        std::unique_ptr<MAST::FieldFunction<RealMatrixX> > _stiff_D;

        std::unique_ptr<MAST::FieldFunction<RealMatrixX> > _damp;

        std::unique_ptr<MAST::FieldFunction<RealMatrixX> > _inertia;

        std::unique_ptr<MAST::FieldFunction<RealMatrixX> > _thermal_A;

        std::unique_ptr<MAST::FieldFunction<RealMatrixX> > _thermal_B;

        std::unique_ptr<MAST::FieldFunction<RealMatrixX> > _transverse_shear;

    };
    
}



#endif // __mast__solid_1d_section_element_property_card__
