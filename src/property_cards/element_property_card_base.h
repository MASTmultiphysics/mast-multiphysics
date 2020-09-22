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

#ifndef __mast__element_property_card_base__
#define __mast__element_property_card_base__


// MAST includes
#include "base/function_set_base.h"
#include "elasticity/bending_operator.h"


namespace MAST
{
    // Forward decleration
    class MaterialPropertyCardBase;
    class ElementBase;
    class GeomElem;
    template <typename ValType> class FieldFunction;
    
    
    enum StrainType {
        LINEAR_STRAIN,
        NONLINEAR_STRAIN
    };
    
    
    
    class ElementPropertyCardBase:
    public MAST::FunctionSetBase {
        
    public:
        ElementPropertyCardBase():
        MAST::FunctionSetBase(),
        _strain_type(MAST::LINEAR_STRAIN),
        _diagonal_mass(false)
        { }
        
        /*!
         *   virtual destructor
         */
        virtual ~ElementPropertyCardBase() { }
        
        /*!
         *   returns the bending model to be used for the element. Should be
         *   reimplemented in the derived classes
         */
        virtual MAST::BendingOperatorType
        bending_model(const MAST::GeomElem& elem) const = 0;
        
        /*!
         *    returns the extra quadrature order (on top of the system) that
         *    this element should use. By default this is zero, and can be
         *    changed by the derived classes
         */
        virtual int extra_quadrature_order(const MAST::GeomElem& elem) const = 0;
        

        virtual std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
        stiffness_A_matrix(const MAST::ElementBase& e) const = 0;

        virtual std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
        stiffness_B_matrix(const MAST::ElementBase& e) const = 0;

        virtual std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
        stiffness_D_matrix(const MAST::ElementBase& e) const = 0;

        virtual std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
        damping_matrix(const MAST::ElementBase& e) const = 0;

        virtual std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
        inertia_matrix(const MAST::ElementBase& e) const = 0;

        virtual std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
        thermal_expansion_A_matrix(const MAST::ElementBase& e) const = 0;

        virtual std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
        thermal_expansion_B_matrix(const MAST::ElementBase& e) const = 0;

        virtual std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
        transverse_shear_stiffness_matrix(const MAST::ElementBase& e) const = 0;

        virtual std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
        prestress_A_matrix( MAST::ElementBase& e) const = 0;

        virtual std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
        prestress_B_matrix( MAST::ElementBase& e) const = 0;
        
        virtual std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
        thermal_conductance_matrix(const MAST::ElementBase& e) const = 0;

        virtual std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
        thermal_capacitance_matrix(const MAST::ElementBase& e) const= 0;

        virtual const MAST::FieldFunction<Real>*
        section(const MAST::ElementBase& e) const = 0;

        
        /*!
         *   return true if the property is isotropic
         */
        virtual bool if_isotropic() const = 0;
        
        
        /*!
         *   return the material property. This needs to be reimplemented
         *   for individual card type, and should be used only for isotropic
         *   cards.
         */
        virtual const MAST::MaterialPropertyCardBase& get_material() const {
            libmesh_error_msg("Must be implemented in derived class; " << __PRETTY_FUNCTION__ << " in " << __FILE__ << " at line number " << __LINE__);
        }
        
        
        /*!
         *   return the material property. This needs to be reimplemented
         *   for individual card type, and should be used only for isotropic
         *   cards.
         * 
         *   Added by DJN. Reference cplusplus.com/forum/beginner/10639
         */
        virtual const MAST::MaterialPropertyCardBase& set_material(MAST::MaterialPropertyCardBase& mat) const {
            libmesh_error_msg("Must be implemented in derived class; " << __PRETTY_FUNCTION__ << " in " << __FILE__ << " at line number " << __LINE__);
        }
        
        
        /*!
         *   vector in the x-y plane of the element. This should not be the same
         *   as the element x-axis.
         *   Only used by 1D sections. Added for polymorphism enhancement.
         *   Added by DJN.
         */
        virtual RealVectorX& y_vector() {
            libmesh_error_msg("Must be implemented in derived class; " << __PRETTY_FUNCTION__ << " in " << __FILE__ << " at line number " << __LINE__);
        }
        

        /*!
         *   constant reference to vector in the x-y plane of the element. 
         *   This should not be the same as the element x-axis.
         *   Only used by 1D sections. Added for polymorphism enhancement.
         *   Added by DJN. 
         */
        virtual const RealVectorX& y_vector() const {
            libmesh_error_msg("Must be implemented in derived class; " << __PRETTY_FUNCTION__ << " in " << __FILE__ << " at line number " << __LINE__);
        }
        
        
        /*!
         * Returns location of shear center of the section
         */
        virtual const libMesh::Point get_shear_center(const libMesh::Point& p, const Real t) const
        {
            libmesh_error_msg("Must be implemented in derived class; " << __PRETTY_FUNCTION__ << " in " << __FILE__ << " at line number " << __LINE__);
        }
        
        
        /*!
         * Returns derivative of location of shear center of the section
         */
        virtual const libMesh::Point get_shear_center_derivative(MAST::FunctionBase& f, const libMesh::Point& p, const Real t)
        {
            libmesh_error_msg("Must be implemented in derived class; " << __PRETTY_FUNCTION__ << " in " << __FILE__ << " at line number " << __LINE__);
        }
        
        
        /*!
         * Returns the location of the centroid of the section
         */
        virtual const libMesh::Point get_centroid(const libMesh::Point& p, const Real t) const
        {
            libmesh_error_msg("Must be implemented in derived class; " << __PRETTY_FUNCTION__ << " in " << __FILE__ << " at line number " << __LINE__);
        }
        
        
        /*!
         * Returns the derivative of the location of the centroid of the section
         */
        virtual const libMesh::Point get_centroid_derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t) const
        {
            libmesh_error_msg("Must be implemented in derived class; " << __PRETTY_FUNCTION__ << " in " << __FILE__ << " at line number " << __LINE__);
        }
        
        
        /*!
         * Returns the points which define the geometry of the cross section.
         */
        virtual const std::vector<libMesh::Point> get_geom_points(const libMesh::Point& p, const Real t, const uint n=201) const
        {
            libmesh_error_msg("Must be implemented in derived class; " << __PRETTY_FUNCTION__ << " in " << __FILE__ << " at line number " << __LINE__);
        }
        
        
        /*!
         * Return the derivative of points which define the geometry of the 
         * cross section.
         */
        virtual const std::vector<libMesh::Point> get_geom_points_derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t, const uint n=201) const 
        {
            libmesh_error_msg("Must be implemented in derived class; " << __PRETTY_FUNCTION__ << " in " << __FILE__ << " at line number " << __LINE__);
        }
        
        
        /*!
         * Returns a list containing vectors of points which define the 
         * geometry of each hole in the cross section.
         */
        virtual const std::vector<std::vector<libMesh::Point>> get_holes_points(const libMesh::Point& p, const Real t, const uint n=201) const
        {
            std::vector<std::vector<libMesh::Point>> empty_hole_list;
            return empty_hole_list;
        }
        
        
        /*!
         * Returns a list containing vectors of points which define the 
         * hole geometry sensitivity to a parameter
         */
        virtual const std::vector<std::vector<libMesh::Point>> get_holes_points_derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t, const uint n=201) const
        {
            std::vector<std::vector<libMesh::Point>> empty_hole_list;
            return empty_hole_list;
        }
        
        
        /*!
         * Returns location of stress evaluation points of the section
         */
        virtual const std::vector<libMesh::Point> get_stress_points(const libMesh::Point& p, const Real t, const libMesh::Point ps) const
        {
            libmesh_error_msg("Must be implemented in derived class; " << __PRETTY_FUNCTION__ << " in " << __FILE__ << " at line number " << __LINE__);
        }
        
        
        /*!
         * Returns location of stress evaluation points of the section
         */
        virtual const std::vector<libMesh::Point> get_stress_points_derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t, const libMesh::Point dps) const
        {
            libmesh_error_msg("Must be implemented in derived class; " << __PRETTY_FUNCTION__ << " in " << __FILE__ << " at line number " << __LINE__);
        }
        
        
        /*!
         * Only used by 1D sections. Added for polymorphism enhancement.
         * 
         * Added by DJN.
         */
        virtual void init() {
            libmesh_error_msg("Must be implemented in derived class; " << __PRETTY_FUNCTION__ << " in " << __FILE__ << " at line number " << __LINE__);
        }
        
        
        /*!
         *   return the material property. This needs to be reimplemented
         *   for individual card type, and should be used only for isotropic
         *   cards.
         * 
         *   Added by DJN. Reference cplusplus.com/forum/beginner/10639
         */
        virtual const MAST::MaterialPropertyCardBase& set_material(MAST::MaterialPropertyCardBase& mat) const {
            libmesh_error_msg("Not Implemented, this needs to be reimplemented for individual card type; In " << __PRETTY_FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__);
        }
        
        
        /*!
         *   vector in the x-y plane of the element. This should not be the same
         *   as the element x-axis.
         *   Only used by 1D sections. Added for polymorphism enhancement.
         *   Added by DJN.
         */
        virtual RealVectorX& y_vector() {
            libmesh_error_msg("Not Implemented, this needs to be reimplemented for individual card type; In " << __PRETTY_FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__);
        }
        

        /*!
         *   constant reference to vector in the x-y plane of the element. 
         *   This should not be the same as the element x-axis.
         *   Only used by 1D sections. Added for polymorphism enhancement.
         *   Added by DJN. 
         */
        virtual const RealVectorX& y_vector() const {
            libmesh_error_msg("Not Implemented, this needs to be reimplemented for individual card type; In " << __PRETTY_FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__);
        }
        
        /*!
         * Only used by 1D sections. Added for polymorphism enhancement.
         * 
         * Added by DJN.
         */
        virtual void init() {
            libmesh_error_msg("Not Implemented, this needs to be reimplemented for individual card type; In " << __PRETTY_FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__);
        }
        
        
        /*!
         *   dimension of the element for which this property is defined
         */
        virtual unsigned int dim() const = 0;
        
        
        /*!
         *    sets the type of strain to be used, which is LINEAR_STRAIN by
         *    default
         */
        void set_strain(MAST::StrainType strain) {
            _strain_type = strain;
        }
        
        
        /*!
         *    returns the type of strain to be used for this element
         */
        const MAST::StrainType strain_type() const {
            return _strain_type;
        }
        
        
        /*!
         *   sets the bending model to be used for the 1D element
         *   Added by DJN to increase section polymorphism
         */
        virtual void set_bending_model(MAST::BendingOperatorType b)  {
            libmesh_error_msg("Not implemented, this needs to be reimplemented for individual card type; In " << __PRETTY_FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__);
        }
        
        
        /*!
         *    sets the mass matrix to be diagonal or consistent
         */
        void set_diagonal_mass_matrix(bool m) {
            _diagonal_mass = m;
        }
        
        
        /*!
         *    returns the type of strain to be used for this element
         */
        bool if_diagonal_mass_matrix() const {
            return _diagonal_mass;
        }
        
        
        /*!
         *    @returns true if the element prestress has been specified, false
         *    otherwise
         */
        virtual bool if_prestressed() const {
            return this->contains("prestress");
        }
        
        
        virtual bool get_warping_only() const
        {
            return _warping_only;
        }
        
    protected:
        
        /*!
         *    type of nonlinear strain to be used for analysis
         */
        MAST::StrainType _strain_type;
        
        /*!
         *    flag to use a diagonal mass matrix. By default, this is false
         */
        bool _diagonal_mass;
        
        bool _warping_only = false;
    };
    
}



#endif // __mast__element_property_card_base__
