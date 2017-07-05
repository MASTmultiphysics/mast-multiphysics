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

// MAST includes
#include "examples/structural/nastran_model_analysis/nastran_model_analysis.h"
#include "examples/structural/nastran_model_analysis/pynastran.h"
#include "examples/structural/nastran_model_analysis/mast_interface.h"
#include "base/nonlinear_system.h"
#include "solver/slepc_eigen_solver.h"
#include "elasticity/structural_modal_eigenproblem_assembly.h"

// libMesh includes
#include "libmesh/exodusII_io.h"


extern libMesh::LibMeshInit* __init;


MAST::NastranModelAnalysis::NastranModelAnalysis() {
    
}


MAST::NastranModelAnalysis::~NastranModelAnalysis() {
    
    delete _model;
}


void
MAST::NastranModelAnalysis::init(libMesh::ElemType etype, bool if_nonlin) {
    
    // create the mesh
    _model   = new MAST::Model(*__init);
    
    
    /////////////////////////////////////////////////
    // Initialize python file
    /////////////////////////////////////////////////
    Py_Initialize();
    
    /////////////////////////////////////////////////
    // run python file
    /////////////////////////////////////////////////
    PyInit_pynastran();

    read_and_initialize_model(_model);

    
    /////////////////////////////////////////////////
    // close Python
    /////////////////////////////////////////////////
    Py_Finalize();

    _model->initialize_after_mesh();
    _model->get_eq_sys().print_info();
}


MAST::Parameter*
MAST::NastranModelAnalysis::get_parameter(const std::string& nm) {
    
    MAST::Parameter* p;
    
    return p;
}



const libMesh::NumericVector<Real>&
MAST::NastranModelAnalysis::solve(bool if_write_output) {
    
    libMesh::MeshBase           &mesh       = _model->get_mesh();
    libMesh::EquationSystems    &eq_sys     = _model->get_eq_sys();
    MAST::NonlinearSystem       &sys        = _model->get_system();
    MAST::SystemInitialization  &sys_init   = _model->get_system_init();
    MAST::PhysicsDisciplineBase &discipline = _model->get_discipline();
    
    
    
    sys.eigen_solver->set_position_of_spectrum(libMesh::LARGEST_MAGNITUDE);
    sys.set_exchange_A_and_B(true);
    sys.set_n_requested_eigenvalues(10);

    
    // create the nonlinear assembly object
    MAST::StructuralModalEigenproblemAssembly   assembly;
    sys.initialize_condensed_dofs(discipline);
    
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

    
    libMesh::NumericVector<Real> *v = nullptr;
    return *v;
}


const libMesh::NumericVector<Real>&
MAST::NastranModelAnalysis::sensitivity_solve(MAST::Parameter& p,
                                              bool if_write_output) {
    
    libMesh::NumericVector<Real> *v = nullptr;
    
    return *v;
}
