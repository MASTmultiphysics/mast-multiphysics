/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2017  Manav Bhatia
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

extern "C" {
#include "Python.h"
}

// C++ includes
#include <string>


// MAST includes
#include "examples/structural/nastran_model_analysis/nastran_model_analysis.h"
#include "examples/structural/nastran_model_analysis/pynastran.h"
#include "examples/structural/nastran_model_analysis/mast_interface.h"
#include "base/nonlinear_system.h"
#include "solver/slepc_eigen_solver.h"
#include "elasticity/structural_nonlinear_assembly.h"
#include "elasticity/structural_modal_eigenproblem_assembly.h"

// libMesh includes
#include "libmesh/exodusII_io.h"
#include "libmesh/dof_map.h"
#include "libmesh/numeric_vector.h"


extern libMesh::LibMeshInit* __init;


MAST::NastranModelAnalysis::NastranModelAnalysis() {
    
    /////////////////////////////////////////////////
    // Initialize python file
    /////////////////////////////////////////////////
    Py_Initialize();
    
}


MAST::NastranModelAnalysis::~NastranModelAnalysis() {
    
    delete _model;

    /////////////////////////////////////////////////
    // close Python
    /////////////////////////////////////////////////
    Py_Finalize();
}


void
MAST::NastranModelAnalysis::init(libMesh::ElemType etype, bool if_nonlin) {
    
    // create the mesh
    std::string
    //f_nm = "/Users/manav/Documents/codes/numerical_lib/pynastran/pyNastran/models/plate/plate.bdf";
    f_nm = "/Users/manav/Documents/codes/tmp/plate.bdf";
    _model   = new MAST::Model(*__init, f_nm);
    
    
    /////////////////////////////////////////////////
    // run python file
    /////////////////////////////////////////////////
    PyInit_pynastran();

    read_and_initialize_model(_model, _model->file_path());

    // make sure that all subcases point to the same SPC set. At this point
    // only a single SPC is handled, but this will likely change in the
    // near future
    std::vector<MAST::Model::SubCase>::const_iterator
    subcase_it  = _model->get_subcases().begin(),
    subcase_end = _model->get_subcases().end();
    
    const unsigned int
    spc_set = subcase_it->spc;
    
    for ( ; subcase_it != subcase_end; subcase_it++) {
        
        libmesh_assert_equal_to(spc_set, subcase_it->spc);
    }
    
    // this sets up the variables and system
    _model->initialize_after_mesh();
    _model->get_eq_sys().print_info();
}


MAST::Parameter*
MAST::NastranModelAnalysis::get_parameter(const std::string& nm) {
    
    MAST::Parameter* p;
    
    libmesh_assert(false); // to be implemented
    
    return p;
}



const libMesh::NumericVector<Real>&
MAST::NastranModelAnalysis::solve(bool if_write_output) {
    
    
    bool
    spc_initialized = false;
    
    // iterate over the set of subcases and process each
    for (unsigned int i=0; i<_model->n_subcases(); i++) {
        
        const MAST::Model::SubCase& s = _model->get_subcase(i);
        std::string nm = _model->file_path();
        read_and_initialize_loads_for_subcase(s.sid, _model, nm);
        
        // now initialize the Dirichlet conditions corresponding to the
        // SPCs. As noted above, this is constrained to a single SPC set
        // shared by all subcases. This will change in future.
        if (!spc_initialized) {
            
            _initialize_dof_constraints();
            spc_initialized = true;
        }
        
        _solve_subcase();
    }
    
    
    libMesh::NumericVector<Real> *v = nullptr;
    return *v;
}


const libMesh::NumericVector<Real>&
MAST::NastranModelAnalysis::sensitivity_solve(MAST::Parameter& p,
                                              bool if_write_output) {
    
    libmesh_assert(false); // to be implemented
    
    libMesh::NumericVector<Real> *v = nullptr;
    
    return *v;
}




void
MAST::NastranModelAnalysis::_solve_subcase() {
    
    
    // make sure
    
    switch (_model->get_sol()) {
        case 101:
            _linear_static_solve();
            break;
            
        case 103:
            _normal_modes_solve();
            break;
            
        default:
            libMesh::out
            << "Sol not implemented: " << _model->get_sol() << std::endl;
            libmesh_assert(false);
    }
}



void
MAST::NastranModelAnalysis::_initialize_dof_constraints() {
    
    std::vector<MAST::Model::SPC>
    &spcs = _model->get_spcs();
    
    boost::bimap<unsigned int, libMesh::Node*>
    &node_id_map = _model->get_node_id_map();
    
    libMesh::DofMap
    &dof_map = _model->get_system().get_dof_map();
    
    // iterate over all the SPCs and tell the dof map about the constraints
    std::vector<MAST::Model::SPC>::const_iterator
    spc_it  = spcs.begin(),
    spc_end = spcs.end();
    
    std::vector<libMesh::dof_id_type> di;
    
    for ( ; spc_it != spc_end; spc_it++) {
        
        // get the node pointer from the map
        const libMesh::Node*
        node    = node_id_map.left.find(spc_it->node)->second;
        
        // get dof number from the node
        di.clear();
        dof_map.dof_indices(node, di, spc_it->comp-1);
        libmesh_assert_equal_to(di.size(), 1);
        
        // now prepare the constraint row and add it to the dof map
        libMesh::DofConstraintRow c_row;
        // libMesh this assumes a 1.0 value for di[0]
        dof_map.add_constraint_row(di[0], c_row, spc_it->val, true);
    }
}





void
MAST::NastranModelAnalysis::_linear_static_solve() {
    
    libMesh::MeshBase           &mesh       = _model->get_mesh();
    libMesh::EquationSystems    &eq_sys     = _model->get_eq_sys();
    MAST::NonlinearSystem       &sys        = _model->get_system();
    MAST::SystemInitialization  &sys_init   = _model->get_system_init();
    MAST::PhysicsDisciplineBase &discipline = _model->get_discipline();
    
    
    // create the nonlinear assembly object
    MAST::StructuralNonlinearAssembly   assembly;
    
    assembly.attach_discipline_and_system(discipline, sys_init);
    
    // zero the solution before solving
    sys.solution->zero();
    
    sys.solve();
    
    // evaluate the outputs
    assembly.calculate_outputs(*(sys.solution));
    
    assembly.clear_discipline_and_system();
    
    // write the solution for visualization
    libMesh::ExodusII_IO(mesh).write_equation_systems("output.exo",
                                                      eq_sys);
}




void
MAST::NastranModelAnalysis::_normal_modes_solve() {
    
    libMesh::MeshBase           &mesh       = _model->get_mesh();
    libMesh::EquationSystems    &eq_sys     = _model->get_eq_sys();
    MAST::NonlinearSystem       &sys        = _model->get_system();
    MAST::SystemInitialization  &sys_init   = _model->get_system_init();
    MAST::PhysicsDisciplineBase &discipline = _model->get_discipline();
    
    
    sys.eigen_solver->set_position_of_spectrum(libMesh::LARGEST_MAGNITUDE);
    sys.set_exchange_A_and_B(true);
    sys.set_n_requested_eigenvalues(10);
    

    // set the constrained dofs and initialize the condensed dof
    sys.initialize_condensed_dofs(discipline);

    
    // create the nonlinear assembly object
    MAST::StructuralModalEigenproblemAssembly   assembly;
    
    assembly.attach_discipline_and_system(discipline, sys_init);
    sys.eigenproblem_solve();
    assembly.clear_discipline_and_system();
    
    // Get the number of converged eigen pairs.
    unsigned int
    nconv = std::min(sys.get_n_converged_eigenvalues(),
                     sys.get_n_requested_eigenvalues());
    
    
    libMesh::ExodusII_IO writer(mesh);
    
    for (unsigned int i=0; i<nconv; i++) {
        
        // now write the eigenvalue
        Real
        re = 0.,
        im = 0.;
        sys.get_eigenpair(i, re, im, *sys.solution);
        
        libMesh::out
        << std::setw(35) << std::fixed << std::setprecision(15)
        << re << std::endl;
        
        // We write the file in the ExodusII format.
        writer.write_timestep("modes.exo",
                              eq_sys,
                              i+1, i);
    }
}



