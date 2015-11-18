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
#include "base/boundary_condition_base.h"


MAST::StressStrainOutputBase::Data::Data(const RealVectorX& stress,
                                         const RealVectorX& strain,
                                         const libMesh::Point& qp,
                                         const libMesh::Point& xyz,
                                         Real JxW):
_stress(stress),
_strain(strain),
_qp(qp),
_xyz(xyz),
_JxW(JxW) {

    // make sure that both the stress and strain are for a 3D configuration,
    // which is the default for this data structure
    libmesh_assert_equal_to(stress.size(), 6);
    libmesh_assert_equal_to(strain.size(), 6);
    
}



const libMesh::Point&
MAST::StressStrainOutputBase::Data::
point_location_in_element_coordinate() const {

    return _qp;
}


const RealVectorX&
MAST::StressStrainOutputBase::Data::stress() const {
    
    return _stress;
}



const RealVectorX&
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


Real
MAST::StressStrainOutputBase::Data::quadrature_point_JxW() const {
    
    return _JxW;
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

    RealVectorX
    dp = RealVectorX::Zero(6);
    
    // if p == 0, then the sensitivity returns nan
    // Hennce, we are avoiding this by setting it to zero whenever p = 0.
    if (fabs(p) > 0.)
        dp =
        (((_dstress_dX.row(0) - _dstress_dX.row(1)) * (_stress(0) - _stress(1)) +
          (_dstress_dX.row(1) - _dstress_dX.row(2)) * (_stress(1) - _stress(2)) +
          (_dstress_dX.row(2) - _dstress_dX.row(0)) * (_stress(2) - _stress(0))) +
         6.0 * (_dstress_dX.row(3) * _stress(3)+
                _dstress_dX.row(4) * _stress(4)+
                _dstress_dX.row(5) * _stress(5))) * 0.5 * pow(p, -0.5);
    
    return dp;
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
           pow(_stress(5), 2)),              //     tau_zz^2)
    dp = 0.;
    
    // if p == 0, then the sensitivity returns nan
    // Hennce, we are avoiding this by setting it to zero whenever p = 0.
    if (fabs(p) > 0.)
        dp =
        (((dstress_dp(0) - dstress_dp(1)) * (_stress(0) - _stress(1)) +
          (dstress_dp(1) - dstress_dp(2)) * (_stress(1) - _stress(2)) +
          (dstress_dp(2) - dstress_dp(0)) * (_stress(2) - _stress(0))) +
         6.0 * (dstress_dp(3) * _stress(3)+
                dstress_dp(4) * _stress(4)+
                dstress_dp(5) * _stress(5))) * 0.5 * pow(p, -0.5);
    
    return dp;
}



MAST::StressStrainOutputBase::StressStrainOutputBase():
MAST::OutputFunctionBase(MAST::STRAIN_STRESS_TENSOR),
_vol_loads(NULL) {
    
}




MAST::StressStrainOutputBase::~StressStrainOutputBase() {
    
    this->clear(true);
}




void
MAST::StressStrainOutputBase::clear(bool clear_elem_subset) {
    
    std::map<const libMesh::Elem*, std::vector<MAST::StressStrainOutputBase::Data*> >::iterator
    map_it  =  _stress_data.begin(),
    map_end =  _stress_data.end();
    
    for ( ; map_it != map_end; map_it++) {
        
        // iterate over all the data and delete them
        std::vector<MAST::StressStrainOutputBase::Data*>::iterator
        it   =  map_it->second.begin(),
        end  =  map_it->second.end();
        
        for ( ; it != end; it++)
            delete *it;
        
        map_it->second.clear();
    }
    
    _stress_data.clear();
    
    if (clear_elem_subset) {
        _elem_subset.clear();
        _vol_loads = NULL;
    }
    
}





bool
MAST::StressStrainOutputBase::
evaluate_for_element(const libMesh::Elem& elem) const {

    bool rval = false;
    
    if ((_elem_subset.size() && _elem_subset.count(&elem)) ||
        (_elem_subset.size() == 0))
        rval = true;
    else
        rval = false;
    
    return rval;
}



void
MAST::StressStrainOutputBase::
set_elements_in_domain(const std::set<const libMesh::Elem*>& elems) {
    
    // make sure that the no data exists
    libmesh_assert(_stress_data.size() == 0);
    libmesh_assert(_elem_subset.size() == 0);
    
    _elem_subset = elems;
}



const std::set<const libMesh::Elem*>&
MAST::StressStrainOutputBase::get_elem_subset() const {
    
    return _elem_subset;
}




MAST::StressStrainOutputBase::Data&
MAST::StressStrainOutputBase::
add_stress_strain_at_qp_location(const libMesh::Elem* e,
                                 const libMesh::Point& quadrature_pt,
                                 const libMesh::Point& physical_pt,
                                 const RealVectorX& stress,
                                 const RealVectorX& strain,
                                 Real JxW) {
    
    // if the element subset has been provided, then make sure that this
    // element exists in the subdomain.
    if (_elem_subset.size())
        libmesh_assert(_elem_subset.count(e));
    
    
    MAST::StressStrainOutputBase::Data* d =
    new MAST::StressStrainOutputBase::Data (stress,
                                            strain,
                                            quadrature_pt,
                                            physical_pt,
                                            JxW);

    
    // check if the specified element exists in the map. If not, add it
    std::map<const libMesh::Elem*, std::vector<MAST::StressStrainOutputBase::Data*> >::iterator
    it = _stress_data.find(e);
    
    // if the element does not exist in the map, add it to the map.
    if (it == _stress_data.end())
        it =
        _stress_data.insert(std::pair<const libMesh::Elem*, std::vector<MAST::StressStrainOutputBase::Data*> >
                            (e, std::vector<MAST::StressStrainOutputBase::Data*>())).first;
    
    it->second.push_back(d);
    
    return *d;
}



const std::map<const libMesh::Elem*, std::vector<MAST::StressStrainOutputBase::Data*> >&
MAST::StressStrainOutputBase::get_stress_strain_data() const {
    
    return _stress_data;
}



unsigned int
MAST::StressStrainOutputBase::
n_elem_in_storage() const {
    
    return (unsigned int)_stress_data.size();
}


unsigned int
MAST::StressStrainOutputBase::
n_stress_strain_data_for_elem(const libMesh::Elem* e) const {
    
    unsigned int n = 0;
    
    std::map<const libMesh::Elem*, std::vector<MAST::StressStrainOutputBase::Data*> >::const_iterator
    it = _stress_data.find(e);
    
    if ( it != _stress_data.end())
        n = (unsigned int)it->second.size();
    
    return n;
}



const std::vector<MAST::StressStrainOutputBase::Data*>&
MAST::StressStrainOutputBase::
get_stress_strain_data_for_elem(const libMesh::Elem *e) const {
    
    // check if the specified element exists in the map. If not, add it
    std::map<const libMesh::Elem*, std::vector<MAST::StressStrainOutputBase::Data*> >::const_iterator
    it = _stress_data.find(e);

    // make sure that the specified elem exists in the map
    libmesh_assert(it != _stress_data.end());
    
    return it->second;
}



void
MAST::StressStrainOutputBase::set_volume_loads(MAST::VolumeBCMapType& vol_loads) {
    
    // make sure that no existing objects are related
    libmesh_assert(!_vol_loads);
    
    _vol_loads  =  &vol_loads;
    
}



MAST::BoundaryConditionBase*
MAST::StressStrainOutputBase::get_thermal_load_for_elem(const libMesh::Elem& elem) {

    MAST::BoundaryConditionBase *rval = NULL;
    
    if (_vol_loads) {
        
        std::pair<std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*>::const_iterator,
        std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*>::const_iterator> it;
        
        it =  _vol_loads->equal_range(elem.subdomain_id());

        for ( ; it.first != it.second; it.first++)
            if (it.first->second->type() == MAST::TEMPERATURE) {
                
                // make sure that only one thermal load exists for an element
                libmesh_assert(!rval);
                
                rval = it.first->second;
            }
    }
    
    
    return rval;
}




Real
MAST::StressStrainOutputBase::
von_Mises_p_norm_functional_for_all_elems(const Real p) const {
    
    Real
    max_val  = 0.,
    e_val    = 0.,
    JxW_val  = 0.,
    JxW      = 0.,
    val      = 0.;
    
    // first find the data with the maximum value, to be used for scaling
    
    // iterate over all element data
    std::map<const libMesh::Elem*, std::vector<MAST::StressStrainOutputBase::Data*> >::const_iterator
    map_it   =  _stress_data.begin(),
    map_end  =  _stress_data.end();
    
    
    for ( ; map_it != map_end; map_it++) {
        
        std::vector<MAST::StressStrainOutputBase::Data*>::const_iterator
        vec_it   = map_it->second.begin(),
        vec_end  = map_it->second.end();
        
        for ( ; vec_it != vec_end; vec_it++) {
            
            // ask this data point for the von Mises stress value
            e_val    =   (*vec_it)->von_Mises_stress();
            
            (e_val > max_val) ?  max_val = e_val: 0; // to find the maximum value
        }
    }
    
    
    // now that we have the maximum value, we evaluate the p-norm
    map_it   =  _stress_data.begin();
    
    for ( ; map_it != map_end; map_it++) {
        
        std::vector<MAST::StressStrainOutputBase::Data*>::const_iterator
        vec_it   = map_it->second.begin(),
        vec_end  = map_it->second.end();
        
        for ( ; vec_it != vec_end; vec_it++) {
            
            // ask this data point for the von Mises stress value
            e_val    =   (*vec_it)->von_Mises_stress();
            JxW      =   (*vec_it)->quadrature_point_JxW();
            
            // we do not use absolute value here, since von Mises stress
            // is >= 0.
            val     +=   pow(e_val/max_val, p) * JxW;
            JxW_val +=   JxW;
        }
    }
    
    val   = max_val * pow(val/JxW_val, 1./p);
    
    return val;
}



Real
MAST::StressStrainOutputBase::
von_Mises_p_norm_functional_sensitivity_for_all_elems
(const Real p,
 const MAST::FunctionBase* f) const {
    
    
    Real
    max_val  = 0.,
    de_val   = 0.,
    e_val    = 0.,
    JxW_val  = 0.,
    JxW      = 0.,
    val      = 0.,
    dval     = 0.;
    
    // first find the data with the maximum value, to be used for scaling
    
    // iterate over all element data
    std::map<const libMesh::Elem*, std::vector<MAST::StressStrainOutputBase::Data*> >::const_iterator
    map_it   =  _stress_data.begin(),
    map_end  =  _stress_data.end();
    
    
    for ( ; map_it != map_end; map_it++) {
        
        std::vector<MAST::StressStrainOutputBase::Data*>::const_iterator
        vec_it   = map_it->second.begin(),
        vec_end  = map_it->second.end();
        
        for ( ; vec_it != vec_end; vec_it++) {
            
            // ask this data point for the von Mises stress value
            e_val    =   (*vec_it)->von_Mises_stress();
            
            (e_val > max_val) ?  max_val = e_val: 0; // to find the maximum value
        }
    }
    
    
    // now that we have the maximum value, we evaluate the p-norm
    map_it   =  _stress_data.begin();
    
    
    for ( ; map_it != map_end; map_it++) {
        
        std::vector<MAST::StressStrainOutputBase::Data*>::const_iterator
        vec_it   = map_it->second.begin(),
        vec_end  = map_it->second.end();
        
        for ( ; vec_it != vec_end; vec_it++) {
            
            // ask this data point for the von Mises stress value
            e_val    =   (*vec_it)->von_Mises_stress();
            de_val   =   (*vec_it)->dvon_Mises_stress_dp(f);
            JxW      =   (*vec_it)->quadrature_point_JxW();
            
            // we do not use absolute value here, since von Mises stress
            // is >= 0.
            val     +=   pow(e_val/max_val, p) * JxW;
            dval    +=   p * pow(e_val/max_val, p-1.) * JxW * de_val/max_val;
            JxW_val +=   JxW;
        }
    }
    
    val   = 1./p * max_val / pow(JxW_val, 1./p) * pow(val, 1./p-1.) * dval;
    
    return val;
}




RealVectorX
MAST::StressStrainOutputBase::
von_Mises_p_norm_functional_state_derivartive_for_all_elems(const Real p) const {
    
    
    Real
    max_val  = 0.,
    e_val    = 0.,
    JxW_val  = 0.,
    JxW      = 0.,
    val      = 0.;
    
    unsigned int
    n_dofs   = 0;
    
    
    // first find the data with the maximum value, to be used for scaling
    
    // iterate over all element data
    std::map<const libMesh::Elem*, std::vector<MAST::StressStrainOutputBase::Data*> >::const_iterator
    map_it   =  _stress_data.begin(),
    map_end  =  _stress_data.end();
    
    
    for ( ; map_it != map_end; map_it++) {
        
        std::vector<MAST::StressStrainOutputBase::Data*>::const_iterator
        vec_it   = map_it->second.begin(),
        vec_end  = map_it->second.end();
        
        for ( ; vec_it != vec_end; vec_it++) {
            
            // ask this data point for the von Mises stress value
            e_val    =   (*vec_it)->von_Mises_stress();
            (e_val > max_val) ?  max_val = e_val: 0; // to find the maximum value
            (n_dofs==0)? n_dofs = (unsigned int)(*vec_it)->dvon_Mises_stress_dX().size(): 0;
        }
    }
    
    
    // now create the vector with the correct dimensions
    RealVectorX
    de_val     = RealVectorX::Zero(n_dofs),
    dval       = RealVectorX::Zero(n_dofs);

    // now that we have the maximum value, we evaluate the p-norm
    map_it   =  _stress_data.begin();
    
    
    for ( ; map_it != map_end; map_it++) {
        
        std::vector<MAST::StressStrainOutputBase::Data*>::const_iterator
        vec_it   = map_it->second.begin(),
        vec_end  = map_it->second.end();
        
        for ( ; vec_it != vec_end; vec_it++) {
            
            // ask this data point for the von Mises stress value
            e_val    =   (*vec_it)->von_Mises_stress();
            de_val   =   (*vec_it)->dvon_Mises_stress_dX();
            JxW      =   (*vec_it)->quadrature_point_JxW();
            
            // we do not use absolute value here, since von Mises stress
            // is >= 0.
            val     +=   pow(e_val/max_val, p) * JxW;
            dval    +=   p * pow(e_val/max_val, p-1.) * JxW * de_val/max_val;
            JxW_val +=   JxW;
        }
    }
    
    return 1./p * max_val / pow(JxW_val, 1./p) * pow(val, 1./p-1.) * dval;
}
