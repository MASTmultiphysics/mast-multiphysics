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


// MAST includes
#include "elasticity/stress_output_base.h"
#include "elasticity/structural_element_base.h"
#include "base/assembly_base.h"
#include "base/physics_discipline_base.h"
#include "base/boundary_condition_base.h"
#include "property_cards/element_property_card_1D.h"
#include "mesh/local_elem_fe.h"
#include "level_set/level_set_intersection.h"


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



void
MAST::StressStrainOutputBase::Data::clear_sensitivity_data() {
    
    _stress_sensitivity.clear();
    _strain_sensitivity.clear();
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
MAST::StressStrainOutputBase::Data::set_sensitivity(const MAST::FunctionBase& f,
                                                    const RealVectorX& dstress_df,
                                                    const RealVectorX& dstrain_df) {

    // make sure that both the stress and strain are for a 3D configuration,
    // which is the default for this data structure
    libmesh_assert_equal_to(dstress_df.size(), 6);
    libmesh_assert_equal_to(dstrain_df.size(), 6);
    
    _stress_sensitivity[&f] = dstress_df;
    _strain_sensitivity[&f] = dstrain_df;
}



const RealVectorX&
MAST::StressStrainOutputBase::Data::
get_stress_sensitivity(const MAST::FunctionBase& f) const {
    
    // make sure that the data exists
    std::map<const MAST::FunctionBase*, RealVectorX>::const_iterator
    it = _stress_sensitivity.find(&f);
    
    libmesh_assert(it != _stress_sensitivity.end());
    
    
    return it->second;
}



const RealVectorX&
MAST::StressStrainOutputBase::Data::
get_strain_sensitivity(const MAST::FunctionBase& f) const {
    
    // make sure that the data exists
    std::map<const MAST::FunctionBase*, RealVectorX>::const_iterator
    it = _strain_sensitivity.find(&f);
    
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
    dp = RealVectorX::Zero(_dstress_dX.cols());
    
    // if p == 0, then the sensitivity returns nan
    // Hence, we are avoiding this by setting it to zero whenever p = 0.
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
dvon_Mises_stress_dp(const MAST::FunctionBase& f) const {
    
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
    // Hence, we are avoiding this by setting it to zero whenever p = 0.
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
MAST::OutputAssemblyElemOperations(),
_p_norm                   (2.),
_rho                      (1.),
_sigma0                   (0.),
_exp_arg_lim              (2.),
_primal_data_initialized  (false),
_JxW_val                  (0.),
_sigma_vm_int             (0.),
_sigma_vm_p_norm          (0.),
_if_stress_plot_mode      (false) {
    
}




MAST::StressStrainOutputBase::~StressStrainOutputBase() {
    
    this->clear();
}




void
MAST::StressStrainOutputBase::zero_for_analysis() {
    
    _JxW_val         = 0.;
    _sigma_vm_int    = 0.;
    _sigma_vm_p_norm = 0.;
}



void
MAST::StressStrainOutputBase::zero_for_sensitivity() {

    // nothign to be done here.
}


void
MAST::StressStrainOutputBase::evaluate() {

    // make sure that this has not been initialized ana calculated for all elems
    libmesh_assert(_physics_elem);
    libmesh_assert(!_primal_data_initialized);

    if (this->if_evaluate_for_element(_physics_elem->elem()))
        // ask for the values
        dynamic_cast<MAST::StructuralElementBase*>
        (_physics_elem)->calculate_stress(false,
                                          nullptr,
                                          *this);
}



void
MAST::StressStrainOutputBase::evaluate_sensitivity(const MAST::FunctionBase &f) {
    
    // the primal data should have been calculated
    libmesh_assert(_physics_elem);
    if (!_if_stress_plot_mode)
        libmesh_assert(_primal_data_initialized);

    if (this->if_evaluate_for_element(_physics_elem->elem())) {
        
        // ask for the values
        dynamic_cast<MAST::StructuralElementBase*>
        (_physics_elem)->calculate_stress(false,
                                          &f,
                                          *this);
    }
}



void
MAST::StressStrainOutputBase::
evaluate_topology_sensitivity(const MAST::FunctionBase &f,
                              const MAST::LevelSetIntersection &intersect,
                              const MAST::FieldFunction<RealVectorX> &vel) {
    
    // the primal data should have been calculated
    libmesh_assert(_physics_elem);
    libmesh_assert(f.is_topology_parameter());
    if (!_if_stress_plot_mode)
        libmesh_assert(_primal_data_initialized);
    
    const libMesh::Elem& elem = _physics_elem->elem();

    // sensitivity only exists at the boundary. So, we proceed with calculation
    // only if this element has an intersection in the interior, or with a side.

    if (this->if_evaluate_for_element(elem)    &&
        intersect.if_elem_has_boundary()       &&
        intersect.has_side_on_interface(elem)) {
        
        // ask for the values
        dynamic_cast<MAST::StructuralElementBase*>
        (_physics_elem)->calculate_stress_boundary_velocity(f, *this,
                                                            intersect.get_side_on_interface(elem),
                                                            vel);
    }
}



Real
MAST::StressStrainOutputBase::output_total() {

    libmesh_assert(!_if_stress_plot_mode);
    
    // if this has not been initialized, then we should do so now
    if (!_primal_data_initialized)
        this->von_Mises_p_norm_functional_for_all_elems();
    
    return _sigma_vm_p_norm;
}


Real
MAST::StressStrainOutputBase::output_sensitivity_for_elem(const MAST::FunctionBase& p) {
    
    libmesh_assert(!_if_stress_plot_mode);
    libmesh_assert(_primal_data_initialized);
    
    Real val = 0.;
    
    this->von_Mises_p_norm_functional_sensitivity_for_elem(p, _physics_elem->elem().id(), val);
    return val;
}



Real
MAST::StressStrainOutputBase::output_sensitivity_total(const MAST::FunctionBase& p) {
    
    libmesh_assert(!_if_stress_plot_mode);
    libmesh_assert(_primal_data_initialized);

    Real
    val   = 0.,
    val_b = 0.;
    
    this->von_Mises_p_norm_functional_sensitivity_for_all_elems(p, val);
    if (p.is_topology_parameter())
        this->von_Mises_p_norm_functional_boundary_sensitivity_for_all_elems(p, val_b);
    return val+val_b;
}



void
MAST::StressStrainOutputBase::output_derivative_for_elem(RealVectorX& dq_dX) {

    libmesh_assert(_physics_elem);
    libmesh_assert(!_if_stress_plot_mode);
    libmesh_assert(_primal_data_initialized);
    
    if (this->if_evaluate_for_element(_physics_elem->elem())) {
        
        dq_dX.setZero();
        
        dynamic_cast<MAST::StructuralElementBase*>
        (_physics_elem)->calculate_stress(true,
                                          nullptr,
                                          *this);
        
        this->von_Mises_p_norm_functional_state_derivartive_for_elem(_physics_elem->elem(),
                                                                     dq_dX);
    }
}





void
MAST::StressStrainOutputBase::init(const libMesh::Elem& elem) {
    
    libmesh_assert(!_physics_elem);
    libmesh_assert(_assembly);
    libmesh_assert(_system);
    
    const MAST::ElementPropertyCardBase& p =
    dynamic_cast<const MAST::ElementPropertyCardBase&>(_discipline->get_property_card(elem));
    
    _physics_elem =
    MAST::build_structural_element(*_system, *_assembly, elem, p).release();
}




void
MAST::StressStrainOutputBase::set_local_fe_data(MAST::LocalElemFE& fe,
                                                const libMesh::Elem& e) const {
    
    if (e.dim() == 1) {
        
        const MAST::ElementPropertyCard1D& p_card =
        dynamic_cast<const MAST::ElementPropertyCard1D&>(_discipline->get_property_card(e));
        
        fe.set_1d_y_vector(p_card.y_vector());
    }
}



void
MAST::StressStrainOutputBase::clear() {
    
    std::map<const libMesh::dof_id_type, std::vector<MAST::StressStrainOutputBase::Data*> >::iterator
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

    
    map_it  =  _boundary_stress_data.begin();
    map_end =  _boundary_stress_data.end();
    
    for ( ; map_it != map_end; map_it++) {
        
        // iterate over all the data and delete them
        std::vector<MAST::StressStrainOutputBase::Data*>::iterator
        it   =  map_it->second.begin(),
        end  =  map_it->second.end();
        
        for ( ; it != end; it++)
            delete *it;
        
        map_it->second.clear();
    }
    
    _boundary_stress_data.clear();

    this->clear_elem();
    
    _primal_data_initialized = false;
    _JxW_val                 = 0.;
    _sigma_vm_int            = 0.;
    _sigma_vm_p_norm         = 0.;
    _if_stress_plot_mode     = false;
}




void
MAST::StressStrainOutputBase::clear_sensitivity_data() {
    
    std::map<const libMesh::dof_id_type, std::vector<MAST::StressStrainOutputBase::Data*> >::iterator
    map_it  =  _stress_data.begin(),
    map_end =  _stress_data.end();
    
    for ( ; map_it != map_end; map_it++) {
        
        // iterate over all the data and delete them
        std::vector<MAST::StressStrainOutputBase::Data*>::iterator
        it   =  map_it->second.begin(),
        end  =  map_it->second.end();
        
        for ( ; it != end; it++)
            (*it)->clear_sensitivity_data();
    }

    
    map_it  =  _boundary_stress_data.begin();
    map_end =  _boundary_stress_data.end();
    
    for ( ; map_it != map_end; map_it++) {
        
        // iterate over all the data and delete them
        std::vector<MAST::StressStrainOutputBase::Data*>::iterator
        it   =  map_it->second.begin(),
        end  =  map_it->second.end();
        
        for ( ; it != end; it++)
            delete *it;
    }
    _boundary_stress_data.clear();
}






MAST::StressStrainOutputBase::Data&
MAST::StressStrainOutputBase::
add_stress_strain_at_qp_location(const libMesh::Elem* e,
                                 const unsigned int qp,
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
    std::map<const libMesh::dof_id_type, std::vector<MAST::StressStrainOutputBase::Data*> >::iterator
    it = _stress_data.find(e->id());
    
    // if the element does not exist in the map, add it to the map.
    if (it == _stress_data.end())
        it =
        _stress_data.insert(std::pair<const libMesh::dof_id_type, std::vector<MAST::StressStrainOutputBase::Data*> >
                            (e->id(), std::vector<MAST::StressStrainOutputBase::Data*>())).first;
    else
        // this assumes that the previous qp data is provided and
        // therefore, this qp number should be == size of the vector.
        libmesh_assert_equal_to(qp, it->second.size());
    
    it->second.push_back(d);
    
    return *d;
}




MAST::StressStrainOutputBase::Data&
MAST::StressStrainOutputBase::
add_stress_strain_at_boundary_qp_location(const libMesh::Elem* e,
                                          const unsigned int s,
                                          const unsigned int qp,
                                          const libMesh::Point& quadrature_pt,
                                          const libMesh::Point& physical_pt,
                                          const RealVectorX& stress,
                                          const RealVectorX& strain,
                                          Real JxW_Vn) {
    
    // if the element subset has been provided, then make sure that this
    // element exists in the subdomain.
    if (_elem_subset.size())
        libmesh_assert(_elem_subset.count(e));
    
    
    MAST::StressStrainOutputBase::Data* d =
    new MAST::StressStrainOutputBase::Data (stress,
                                            strain,
                                            quadrature_pt,
                                            physical_pt,
                                            JxW_Vn);
    
    
    // check if the specified element exists in the map. If not, add it
    std::map<const libMesh::dof_id_type, std::vector<MAST::StressStrainOutputBase::Data*> >::iterator
    it = _boundary_stress_data.find(e->id());
    
    // if the element does not exist in the map, add it to the map.
    if (it == _boundary_stress_data.end())
        it =
        _boundary_stress_data.insert(std::pair<const libMesh::dof_id_type, std::vector<MAST::StressStrainOutputBase::Data*> >
                                     (e->id(), std::vector<MAST::StressStrainOutputBase::Data*>())).first;
    else
        // this assumes that the previous qp data is provided and
        // therefore, this qp number should be == size of the vector.
        libmesh_assert_equal_to(qp, it->second.size());
    
    it->second.push_back(d);
    
    return *d;
}




const std::map<const libMesh::dof_id_type, std::vector<MAST::StressStrainOutputBase::Data*> >&
MAST::StressStrainOutputBase::get_stress_strain_data() const {
    
    return _stress_data;
}



unsigned int
MAST::StressStrainOutputBase::
n_stress_strain_data_for_elem(const libMesh::Elem* e) const {
    
    unsigned int n = 0;
    
    std::map<const libMesh::dof_id_type, std::vector<MAST::StressStrainOutputBase::Data*> >::const_iterator
    it = _stress_data.find(e->id());
    
    if ( it != _stress_data.end())
        n = (unsigned int)it->second.size();
    
    return n;
}



unsigned int
MAST::StressStrainOutputBase::
n_boundary_stress_strain_data_for_elem(const libMesh::Elem* e) const {
    
    unsigned int n = 0;
    
    std::map<const libMesh::dof_id_type, std::vector<MAST::StressStrainOutputBase::Data*> >::const_iterator
    it = _boundary_stress_data.find(e->id());
    
    if ( it != _boundary_stress_data.end())
        n = (unsigned int)it->second.size();
    
    return n;
}



const std::vector<MAST::StressStrainOutputBase::Data*>&
MAST::StressStrainOutputBase::
get_stress_strain_data_for_elem(const libMesh::Elem *e) const {
    
    // check if the specified element exists in the map. If not, add it
    std::map<const libMesh::dof_id_type, std::vector<MAST::StressStrainOutputBase::Data*> >::const_iterator
    it = _stress_data.find(e->id());

    // make sure that the specified elem exists in the map
    libmesh_assert(it != _stress_data.end());
    
    return it->second;
}



MAST::StressStrainOutputBase::Data&
MAST::StressStrainOutputBase::
get_stress_strain_data_for_elem_at_qp(const libMesh::Elem* e,
                                      const unsigned int qp) {

    
    // check if the specified element exists in the map. If not, add it
    std::map<const libMesh::dof_id_type, std::vector<MAST::StressStrainOutputBase::Data*> >::const_iterator
    it = _stress_data.find(e->id());
    
    // make sure that the specified elem exists in the map
    libmesh_assert(it != _stress_data.end());
    
    libmesh_assert_less(qp, it->second.size());
    
    return *(it->second[qp]);
}



MAST::BoundaryConditionBase*
MAST::StressStrainOutputBase::get_thermal_load_for_elem(const libMesh::Elem& elem) {

    MAST::BoundaryConditionBase *rval = nullptr;
 
    // this should only be called if the user has specified evaluation of
    // stress for this
    libmesh_assert(this->if_evaluate_for_element(elem));
        

    MAST::VolumeBCMapType&
    vol_loads = _discipline->volume_loads();
    
    std::pair<std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*>::const_iterator,
    std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*>::const_iterator> it;
    
    it =  vol_loads.equal_range(elem.subdomain_id());
    
    // FIXME: that this assumes that only one temperature boundary condition
    // is applied for an element.
    for ( ; it.first != it.second; it.first++)
        if (it.first->second->type() == MAST::TEMPERATURE) {
            
            // make sure that only one thermal load exists for an element
            libmesh_assert(!rval);
            
            rval = it.first->second;
        }
    
    return rval;
}




void
MAST::StressStrainOutputBase::
von_Mises_p_norm_functional_for_all_elems() {
    
    libmesh_assert(!_if_stress_plot_mode);
    libmesh_assert(!_primal_data_initialized);
    libmesh_assert_greater(_sigma0, 0.);
    
    Real
    sp               = 0.,
    exp_sp           = 0.,
    e_val            = 0.,
    JxW              = 0.;
    
    _JxW_val         = 0.;
    _sigma_vm_int    = 0.;
    _sigma_vm_p_norm = 0.;
    
    // first find the data with the maximum value, to be used for scaling
    
    // iterate over all element data
    std::map<const libMesh::dof_id_type, std::vector<MAST::StressStrainOutputBase::Data*> >::const_iterator
    map_it   =  _stress_data.begin(),
    map_end  =  _stress_data.end();
    
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
            sp              =  pow(e_val/_sigma0, _p_norm);
            if (_rho * sp > _exp_arg_lim)
                exp_sp          =  exp(_exp_arg_lim);
            else
                exp_sp          =  exp(_rho * sp);
            _sigma_vm_int  +=  sp * exp_sp * JxW;
            _JxW_val       +=  exp_sp * JxW;
        }
    }
    
    
    _sigma_vm_p_norm         = _sigma0 * pow(_sigma_vm_int/_JxW_val, 1./_p_norm);
    _primal_data_initialized = true;
}



void
MAST::StressStrainOutputBase::
von_Mises_p_norm_functional_sensitivity_for_all_elems
(const MAST::FunctionBase& f,
 Real& dsigma_vm_val_df) const {
    
    libmesh_assert(!_if_stress_plot_mode);
    libmesh_assert(_primal_data_initialized);

    Real
    val      = 0.;
    
    dsigma_vm_val_df = 0.;
    
    // iterate over all element data
    std::map<const libMesh::dof_id_type, std::vector<MAST::StressStrainOutputBase::Data*> >::const_iterator
    map_it   =  _stress_data.begin(),
    map_end  =  _stress_data.end();
    
    for ( ; map_it != map_end; map_it++) {
        
        this->von_Mises_p_norm_functional_sensitivity_for_elem(f, map_it->first, val);
        dsigma_vm_val_df += val;
    }
}




void
MAST::StressStrainOutputBase::
von_Mises_p_norm_functional_boundary_sensitivity_for_all_elems
(const MAST::FunctionBase& f,
 Real& dsigma_vm_val_df) const {
    
    libmesh_assert(!_if_stress_plot_mode);
    libmesh_assert(_primal_data_initialized);
    
    Real
    val      = 0.;
    
    dsigma_vm_val_df = 0.;
    
    // iterate over all element data
    std::map<const libMesh::dof_id_type, std::vector<MAST::StressStrainOutputBase::Data*> >::const_iterator
    map_it   =  _boundary_stress_data.begin(),
    map_end  =  _boundary_stress_data.end();
    
    for ( ; map_it != map_end; map_it++) {
        
        this->von_Mises_p_norm_functional_boundary_sensitivity_for_elem(f, map_it->first, val);
        dsigma_vm_val_df += val;
    }
}



void
MAST::StressStrainOutputBase::
von_Mises_p_norm_functional_sensitivity_for_elem
(const MAST::FunctionBase& f,
 const libMesh::dof_id_type e_id,
 Real& dsigma_vm_val_df) const {
    
    libmesh_assert(!_if_stress_plot_mode);
    libmesh_assert(_primal_data_initialized);
    libmesh_assert_greater(_sigma0, 0.);
    
    Real
    sp            = 0.,
    sp1           = 0.,
    exp_sp        = 0.,
    e_val         = 0.,
    de_val        = 0.,
    num_sens      = 0.,
    denom_sens    = 0.,
    JxW           = 0.;
    
    dsigma_vm_val_df = 0.;
    
    // iterate over all element data
    std::map<const libMesh::dof_id_type, std::vector<MAST::StressStrainOutputBase::Data*> >::const_iterator
    map_it   =  _stress_data.find(e_id),
    map_end  =  _stress_data.end();
    
    libmesh_assert(map_it != map_end);
    
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
        sp           =  pow(e_val/_sigma0, _p_norm);
        sp1          =  pow(e_val/_sigma0, _p_norm-1.);
        if (_rho*sp > _exp_arg_lim) {
            exp_sp       =  exp(_exp_arg_lim);
            denom_sens  +=  0.;
            num_sens    += _p_norm * sp1 * exp_sp * de_val/_sigma0 * JxW;
        }
        else {
            exp_sp       =  exp(_rho * sp);
            denom_sens  +=  exp_sp * _rho * _p_norm * sp1 * de_val/_sigma0 * JxW;
            num_sens    += (1. + sp * _rho) * _p_norm * sp1 * exp_sp * de_val/_sigma0 * JxW;
        }
            
    }
    
    dsigma_vm_val_df = _sigma0/_p_norm * pow(_sigma_vm_int/_JxW_val, 1./_p_norm - 1.) *
    (num_sens / _JxW_val - _sigma_vm_int / pow(_JxW_val, 2.) * denom_sens);
}





void
MAST::StressStrainOutputBase::
von_Mises_p_norm_functional_boundary_sensitivity_for_elem
(const MAST::FunctionBase& f,
 const libMesh::dof_id_type e_id,
 Real& dsigma_vm_val_df) const {
    
    libmesh_assert(!_if_stress_plot_mode);
    libmesh_assert(_primal_data_initialized);
    libmesh_assert_greater(_sigma0, 0.);
    
    Real
    sp            = 0.,
    exp_sp        = 0.,
    e_val         = 0.,
    JxW_Vn        = 0.,
    num_sens      = 0.,
    denom_sens    = 0.;

    dsigma_vm_val_df = 0.;
    
    // iterate over all element data
    std::map<const libMesh::dof_id_type, std::vector<MAST::StressStrainOutputBase::Data*> >::const_iterator
    map_it   =  _boundary_stress_data.find(e_id),
    map_end  =  _boundary_stress_data.end();
    
    libmesh_assert(map_it != map_end);
    
    std::vector<MAST::StressStrainOutputBase::Data*>::const_iterator
    vec_it   = map_it->second.begin(),
    vec_end  = map_it->second.end();
    
    for ( ; vec_it != vec_end; vec_it++) {
        
        // ask this data point for the von Mises stress value
        e_val    =   (*vec_it)->von_Mises_stress();
        JxW_Vn   =   (*vec_it)->quadrature_point_JxW();
        
        // we do not use absolute value here, since von Mises stress
        // is >= 0.
        sp           =  pow(e_val/_sigma0, _p_norm);
        if (_rho*sp > _exp_arg_lim)
            exp_sp       =  exp(_exp_arg_lim);
        else
            exp_sp       =  exp(_rho * sp);
        denom_sens  +=  exp_sp * JxW_Vn;
        num_sens    +=  sp * exp_sp * JxW_Vn;
    }
    
    dsigma_vm_val_df = _sigma0/_p_norm * pow(_sigma_vm_int/_JxW_val, 1./_p_norm - 1.) *
    (num_sens / _JxW_val - _sigma_vm_int / pow(_JxW_val, 2.) * denom_sens);
}




void
MAST::StressStrainOutputBase::
von_Mises_p_norm_functional_state_derivartive_for_elem(const libMesh::Elem& e,
                                                       RealVectorX& dq_dX) const {
    
    libmesh_assert(!_if_stress_plot_mode);
    libmesh_assert(_primal_data_initialized);
    libmesh_assert_greater(_sigma0, 0.);
    
    Real
    sp            = 0.,
    sp1           = 0.,
    exp_sp        = 0.,
    e_val         = 0.,
    JxW           = 0.;
    
    RealVectorX
    denom_sens     = RealVectorX::Zero(dq_dX.size()),
    num_sens       = RealVectorX::Zero(dq_dX.size()),
    de_val         = RealVectorX::Zero(dq_dX.size());
    dq_dX.setZero();
    
    // first find the data with the maximum value, to be used for scaling
    
    // iterate over all element data
    std::map<const libMesh::dof_id_type, std::vector<MAST::StressStrainOutputBase::Data*> >::const_iterator
    map_it   =  _stress_data.find(e.id()),
    map_end  =  _stress_data.end();

    // make sure that the data exists
    libmesh_assert(map_it != map_end);
    
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
        sp           =  pow(e_val/_sigma0, _p_norm);
        sp1          =  pow(e_val/_sigma0, _p_norm-1.);
        if (_rho*sp > _exp_arg_lim) {
            
            exp_sp       =  exp(_exp_arg_lim);
            denom_sens  +=  0. * de_val;
            num_sens    +=  1. * _p_norm * sp1 * exp_sp/_sigma0 * JxW * de_val;
        }
        else {
            
            exp_sp       =  exp(_rho * sp);
            denom_sens  +=  exp_sp * _rho * _p_norm * sp1/_sigma0 * JxW * de_val;
            num_sens    += (1. + sp * _rho) * _p_norm * sp1 * exp_sp/_sigma0 * JxW * de_val;
        }
    }
    
    dq_dX = _sigma0/_p_norm * pow(_sigma_vm_int/_JxW_val, 1./_p_norm - 1.) *
    (num_sens / _JxW_val - _sigma_vm_int / pow(_JxW_val, 2.) * denom_sens);
}




