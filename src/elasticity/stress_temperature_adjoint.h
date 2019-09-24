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

#ifndef __mast__stress_temperature_adjoint_h__
#define __mast__stress_temperature_adjoint_h__


// MAST includes
#include "elasticity/stress_output_base.h"


// libMesh includes
#include "libmesh/elem.h"

namespace MAST {
    
    
    /*!
     *  The stress and thermoelastic analysis are dependent on temperature.
     *  This class provides the contribution to the forcing function for
     *  calculation of adjoint vector for the thermal system, as a result
     *  of the stress and themroelstic analyses.
     */
    class StressTemperatureAdjoint:
    public MAST::StressStrainOutputBase {
        
    public:
        
        /*!
         *    default constructor
         */
        StressTemperatureAdjoint(MAST::StressStrainOutputBase& stress);
        
        virtual ~StressTemperatureAdjoint();

        void
        set_thermal_assembly(MAST::AssemblyBase& thermal_assembly);
        
        void
        set_structural_solutions(const libMesh::NumericVector<Real>& sol,
                                 const libMesh::NumericVector<Real>& adj_sol);
        
        /*!
         *    checks to see if the object has been told about the subset of
         *    elements and if the specified element is in the subset.
         */
        virtual bool if_evaluate_for_element(const MAST::GeomElem& elem) const {
            return _stress.if_evaluate_for_element(elem);
        }
        

        
        /*!
         *   sets the element solution
         */
        virtual void
        set_elem_solution(const RealVectorX& sol);

        virtual void output_derivative_for_elem(RealVectorX& dq_dX);
        

        virtual MAST::StressStrainOutputBase::Data&
        add_stress_strain_at_qp_location(const MAST::GeomElem& e,
                                         const unsigned int qp,
                                         const libMesh::Point& quadrature_pt,
                                         const libMesh::Point& physical_pt,
                                         const RealVectorX& stress,
                                         const RealVectorX& strain,
                                         Real JxW) {
            libmesh_error(); // shoudl not get called
        }
        

        /*!
         *   add the stress tensor associated with the \p qp on side \p s of
         *   element \p e. @returns a reference to \p Data.
         */
        virtual MAST::StressStrainOutputBase::Data&
        add_stress_strain_at_boundary_qp_location(const MAST::GeomElem& e,
                                                  const unsigned int s,
                                                  const unsigned int qp,
                                                  const libMesh::Point& quadrature_pt,
                                                  const libMesh::Point& physical_pt,
                                                  const RealVectorX& stress,
                                                  const RealVectorX& strain,
                                                  Real JxW_Vn) {
            libmesh_error(); // should not get called
        }

        /*!
         *    @returns the vector of stress/strain data for specified elem at
         *    the specified quadrature point.
         */
        virtual MAST::StressStrainOutputBase::Data&
        get_stress_strain_data_for_elem_at_qp(const MAST::GeomElem& e,
                                              const unsigned int qp) {
            libmesh_error(); // should not get called
        }
        
        /*!
         *    @returns the map of stress/strain data for all elems
         */
        virtual const std::map<const libMesh::dof_id_type,
        std::vector<MAST::StressStrainOutputBase::Data*> >&
        get_stress_strain_data() const {
            libmesh_error(); // should not get called
        }
        
        
        /*!
         *    @returns the vector of stress/strain data for specified elem.
         */
        virtual const std::vector<MAST::StressStrainOutputBase::Data*>&
        get_stress_strain_data_for_elem(const MAST::GeomElem& e) const {
            libmesh_error(); // should not get called
        }

    protected:

        MAST::StressStrainOutputBase&                  _stress;
        MAST::AssemblyBase*                            _thermal_assembly;
        std::unique_ptr<libMesh::NumericVector<Real>>  _structural_sol;
        std::unique_ptr<libMesh::NumericVector<Real>>  _structural_adjoint;
    };
}

#endif // __mast__stress_temperature_adjoint_h__

