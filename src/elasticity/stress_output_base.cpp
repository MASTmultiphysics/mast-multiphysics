/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2015  Manav Bhatia
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


// MAST includes
#include "elasticity/stress_output_base.h"



MAST::StressStrainOutputBase::Data::Data(const RealVectorX& stress,
                                         const RealVectorX& strain,
                                         const libMesh::Point& qp,
                                         const libMesh::Point& xyz):
_stress(stress),
_strain(strain),
_qp(qp),
_xyz(xyz) {

    // make sure that both the stress and strain are for a 3D configuration,
    // which is the default for this data structure
    libmesh_assert_equal_to(stress.size(), 6);
    libmesh_assert_equal_to(strain.size(), 6);
    
}




const RealVectorX
MAST::StressStrainOutputBase::Data::stress() const {
    
    return _stress;
}



const RealVectorX
MAST::StressStrainOutputBase::Data::strain() const {
    return _strain;
}



void
MAST::StressStrainOutputBase::Data::set_derivatives(const RealMatrixX& dstress_dX,
                                                    const RealMatrixX& dstrain_dX) {
    
    // make sure that the number of rows is 6.
    libmesh_assert_equal_to(dstress_dX.rows(), 6);
    libmesh_assert_equal_to(dstrain_dX.rows(), 6);
    
    _dstress_dX = dstress_dX;
    _dstrain_dX = dstrain_dX;
}



const RealMatrixX&
MAST::StressStrainOutputBase::Data::get_dstress_dX() const {
    
    return _dstress_dX;
}


const RealMatrixX&
MAST::StressStrainOutputBase::Data::get_dstrain_dX() const {
    
    return _dstrain_dX;
}



void
MAST::StressStrainOutputBase::Data::set_sensitivity(const MAST::FunctionBase* f,
                                                    const RealVectorX& dstress_df,
                                                    const RealVectorX& dstrain_df) {

    // make sure that both the stress and strain are for a 3D configuration,
    // which is the default for this data structure
    libmesh_assert_equal_to(dstress_df.size(), 6);
    libmesh_assert_equal_to(dstrain_df.size(), 6);
    
    _stress_sensitivity[f] = dstress_df;
    _strain_sensitivity[f] = dstrain_df;
}



const RealVectorX&
MAST::StressStrainOutputBase::Data::
get_stress_sensitivity(const MAST::FunctionBase* f) const {
    
    // make sure that the data exists
    std::map<const MAST::FunctionBase*, RealVectorX>::const_iterator
    it = _stress_sensitivity.find(f);
    
    libmesh_assert(it != _stress_sensitivity.end());
    
    
    return it->second;
}



const RealVectorX&
MAST::StressStrainOutputBase::Data::
get_strain_sensitivity(const MAST::FunctionBase* f) const {
    
    // make sure that the data exists
    std::map<const MAST::FunctionBase*, RealVectorX>::const_iterator
    it = _strain_sensitivity.find(f);
    
    libmesh_assert(it != _strain_sensitivity.end());
    
    
    return it->second;
}



Real
MAST::StressStrainOutputBase::Data::von_Mises_stress() const {
    
    // make sure that the data is available
    libmesh_assert_equal_to(_stress.size(), 6);

    
    return
    pow(0.5 * (pow(_stress(0)-_stress(1),2) +    //(((sigma_xx - sigma_yy)^2    +
               pow(_stress(1)-_stress(2),2) +    //  (sigma_yy - sigma_zz)^2    +
               pow(_stress(2)-_stress(0),2)) +   //  (sigma_zz - sigma_xx)^2)/2 +
        3.0 * (pow(_stress(3), 2) +              // 3* (tau_xx^2 +
               pow(_stress(4), 2) +              //     tau_yy^2 +
               pow(_stress(5), 2)), 0.5);        //     tau_zz^2))^.5
}



RealVectorX
MAST::StressStrainOutputBase::Data::dvon_Mises_stress_dX() const {
    
    // make sure that the data is available
    libmesh_assert_equal_to(_stress.size(), 6);
    libmesh_assert_equal_to(_dstress_dX.rows(), 6);
    
    Real
    p =
    0.5 * (pow(_stress(0)-_stress(1),2) +    //((sigma_xx - sigma_yy)^2    +
           pow(_stress(1)-_stress(2),2) +    // (sigma_yy - sigma_zz)^2    +
           pow(_stress(2)-_stress(0),2)) +   // (sigma_zz - sigma_xx)^2)/2 +
    3.0 * (pow(_stress(3), 2) +              // 3* (tau_xx^2 +
           pow(_stress(4), 2) +              //     tau_yy^2 +
           pow(_stress(5), 2));              //     tau_zz^2)
    
    return
    (((_dstress_dX.row(0) - _dstress_dX.row(1)) * (_stress(0) - _stress(1)) +
      (_dstress_dX.row(1) - _dstress_dX.row(2)) * (_stress(1) - _stress(2)) +
      (_dstress_dX.row(2) - _dstress_dX.row(0)) * (_stress(2) - _stress(0))) +
     6.0 * (_dstress_dX.row(3) * _stress(3)+
            _dstress_dX.row(4) * _stress(4)+
            _dstress_dX.row(5) * _stress(5))) * 0.5 * pow(p, -0.5);
    
}




Real
MAST::StressStrainOutputBase::Data::
dvon_Mises_stress_dp(const MAST::FunctionBase* f) const {
    
    // make sure that the data is available
    libmesh_assert_equal_to(_stress.size(), 6);
    
    // get the stress sensitivity data
    const RealVectorX dstress_dp = this->get_stress_sensitivity(f);
    
    
    Real
    p =
    0.5 * (pow(_stress(0)-_stress(1),2) +    //((sigma_xx - sigma_yy)^2    +
           pow(_stress(1)-_stress(2),2) +    // (sigma_yy - sigma_zz)^2    +
           pow(_stress(2)-_stress(0),2)) +   // (sigma_zz - sigma_xx)^2)/2 +
    3.0 * (pow(_stress(3), 2) +              // 3* (tau_xx^2 +
           pow(_stress(4), 2) +              //     tau_yy^2 +
           pow(_stress(5), 2));              //     tau_zz^2)

    
    
    return
    (((dstress_dp(0) - dstress_dp(1)) * (_stress(0) - _stress(1)) +
      (dstress_dp(1) - dstress_dp(2)) * (_stress(1) - _stress(2)) +
      (dstress_dp(2) - dstress_dp(0)) * (_stress(2) - _stress(0))) +
     6.0 * (dstress_dp(3) * _stress(3)+
            dstress_dp(4) * _stress(4)+
            dstress_dp(5) * _stress(5))) * 0.5 * pow(p, -0.5);
    
}



MAST::StressStrainOutputBase::StressStrainOutputBase():
MAST::OutputFunctionBase(MAST::STRAIN_STRESS_TENSOR) {
    
}




MAST::StressStrainOutputBase::~StressStrainOutputBase() {
    
    this->clear();
}




void
MAST::StressStrainOutputBase::clear() {
    
    // iterate over all the data and delete them
    std::vector<MAST::StressStrainOutputBase::Data*>::iterator
    it   =  _stress_data.begin(),
    end  =  _stress_data.end();
    
    for ( ; it != end; it++) {
        
        delete *it;
    }
    
    _stress_data.clear();
}





MAST::StressStrainOutputBase::Data& 
MAST::StressStrainOutputBase::
add_stress_strain_at_qp_location(const libMesh::Point& quadrature_pt,
                                 const libMesh::Point& physical_pt,
                                 const RealVectorX& stress,
                                 const RealVectorX& strain) {
    
    
    MAST::StressStrainOutputBase::Data* d =
    new MAST::StressStrainOutputBase::Data (stress,
                                            strain,
                                            quadrature_pt,
                                            physical_pt);
    
    _stress_data.push_back(d);
    
    return *d;
}



const std::vector<MAST::StressStrainOutputBase::Data*>&
MAST::StressStrainOutputBase::get_stress_strain_data() const {
    
    return _stress_data;
}
