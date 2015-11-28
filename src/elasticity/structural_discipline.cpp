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
#include "elasticity/structural_discipline.h"
#include "elasticity/stress_output_base.h"
#include "base/parameter.h"

// libMesh includes
#include "libmesh/explicit_system.h"
#include "libmesh/fe_type.h"
#include "libmesh/mesh_base.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/exodusII_io.h"


MAST::StructuralDiscipline::
StructuralDiscipline(libMesh::EquationSystems& eq_sys):
MAST::PhysicsDisciplineBase(eq_sys)
{ }



MAST::StructuralDiscipline::~StructuralDiscipline()
{ }


void
get_max_stress_strain_values(const std::vector<MAST::StressStrainOutputBase::Data*>& data,
                             RealVectorX&           max_strain,
                             RealVectorX&           max_stress,
                             Real&                  max_vm,
                             const MAST::Parameter* p) {
    
    max_strain    = RealVectorX::Zero(6);
    max_stress    = RealVectorX::Zero(6);
    max_vm        = 0.;
    
    // if there is only one data point, the simply copy the value to the output
    // routines
    if (data.size() == 1) {
        if (p == NULL) {
            max_strain  = data[0]->strain();
            max_stress  = data[0]->stress();
            max_vm      = data[0]->von_Mises_stress();
        }
        else {
            max_strain  = data[0]->get_strain_sensitivity(p);
            max_stress  = data[0]->get_stress_sensitivity(p);
            max_vm      = data[0]->dvon_Mises_stress_dp  (p);
        }
        
        return;
    }
    
    // if multiple values are provided for an element, then we need to compare
    std::vector<MAST::StressStrainOutputBase::Data*>::const_iterator
    it        = data.begin(),
    end       = data.end();
    
    Real
    vm        = 0.;
    
    for ( ; it != end; it++) {
        
        // get the strain value at this point
        const RealVectorX& strain =  (*it)->strain();
        const RealVectorX& stress =  (*it)->stress();
        vm                        =  (*it)->von_Mises_stress();
        
        // now compare
        if (vm > max_vm)                      max_vm        = vm;
        
        for ( unsigned int i=0; i<6; i++) {
            if (fabs(strain(i)) > max_strain(i))  max_strain(i) = strain(i);
            if (fabs(stress(i)) > max_stress(i))  max_stress(i) = stress(i);
        }
    }
}


template <typename ValType>
void MAST::StructuralDiscipline::
plot_stress_strain_data(const std::string& file_nm,
                        const MAST::Parameter* p) const {
    
    libMesh::MeshBase& mesh =  _eq_systems.get_mesh();
    
    // create a new equation system with an explicit system that has one
    // variable, and output the stress data
    libMesh::EquationSystems eq_sys(mesh);
    
    libMesh::ExplicitSystem&
    sys       = eq_sys.add_system<libMesh::ExplicitSystem>("StressOutput");
    
    libMesh::FEType
    fetype(libMesh::CONSTANT, libMesh::MONOMIAL); // constant value per element
    
    std::vector<unsigned int> vars(13);
    
    vars[0]  = sys.add_variable("epsilon-xx", fetype);
    vars[1]  = sys.add_variable("epsilon-yy", fetype);
    vars[2]  = sys.add_variable("epsilon-zz", fetype);
    vars[3]  = sys.add_variable("epsilon-xy", fetype);
    vars[4]  = sys.add_variable("epsilon-yz", fetype);
    vars[5]  = sys.add_variable("epsilon-zx", fetype);
    
    vars[6]  = sys.add_variable("sigma-xx",    fetype);
    vars[7]  = sys.add_variable("sigma-yy",    fetype);
    vars[8]  = sys.add_variable("sigma-zz",    fetype);
    vars[9]  = sys.add_variable("sigma-xy",    fetype);
    vars[10] = sys.add_variable("sigma-yz",    fetype);
    vars[11] = sys.add_variable("sigma-zx",    fetype);

    vars[12] = sys.add_variable("sigma-vm",    fetype);

    eq_sys.init();

    // iterate over the data in the output map and write it
    MAST::VolumeOutputMapType::const_iterator
    it   =  _vol_output_map.begin(),
    end  =  _vol_output_map.end();

    
    libMesh::subdomain_id_type
    sid      = 0;
    
    libMesh::dof_id_type
    dof_id   = 0;

    RealVectorX
    max_strain_vals  =   RealVectorX::Zero(6),
    max_stress_vals  =   RealVectorX::Zero(6);

    Real
    max_vm_stress    =   0.;
    
    
    for ( ; it != end; it++) {
    
        sid      = it->first;
        
        // get reference to the output object
        const MAST::StressStrainOutputBase&
        output   =  dynamic_cast<MAST::StressStrainOutputBase&>(*(it->second));
        
        // get the stress-strain data map from the object
        const std::map<const libMesh::Elem*,
        std::vector<MAST::StressStrainOutputBase::Data*> >& output_map =
        output.get_stress_strain_data();
        
        // now iteragtove over all the elements and set the value in the
        // new system used for output
        std::map<const libMesh::Elem*,
        std::vector<MAST::StressStrainOutputBase::Data*> >::const_iterator
        e_it    =  output_map.begin(),
        e_end   =  output_map.end();
        
        for ( ; e_it != e_end; e_it++) {
            
            get_max_stress_strain_values(e_it->second,
                                         max_strain_vals,
                                         max_stress_vals,
                                         max_vm_stress,
                                         p);
            
            // set the values in the system
            // stress value
            dof_id     =   (e_it->first)->dof_number(sys.number(), vars[12], 0);
            sys.solution->set(dof_id, max_vm_stress);
            
            for (unsigned int i=0; i<6; i++) {
                // strain value
                dof_id     =   (e_it->first)->dof_number(sys.number(), vars[i], 0);
                sys.solution->set(dof_id, max_strain_vals(i));
                
                // stress value
                dof_id     =   (e_it->first)->dof_number(sys.number(), vars[i+6], 0);
                sys.solution->set(dof_id, max_stress_vals(i));
            }
        }
    }
    
    sys.solution->close();
    
    // now output
    ValType(mesh).write_equation_systems(file_nm, eq_sys);
}



// explicit instantiation
template void MAST::StructuralDiscipline::
plot_stress_strain_data<libMesh::ExodusII_IO>(const std::string&     file_nm,
                                              const MAST::Parameter* p) const;

