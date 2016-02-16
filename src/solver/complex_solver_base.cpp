/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2016  Manav Bhatia
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
#include "solver/complex_solver_base.h"
#include "base/complex_assembly_base.h"


// libMesh includes
#include "libmesh/numeric_vector.h"
#include "libmesh/system.h"



MAST::ComplexSolverBase::ComplexSolverBase():
_assembly(NULL),
tol(1.0e-3),
max_iters(20) {
    
}



MAST::ComplexSolverBase::~ComplexSolverBase() {
    
}




void
MAST::ComplexSolverBase::set_assembly(MAST::ComplexAssemblyBase& assembly) {
    
    _assembly = &assembly;
}



void
MAST::ComplexSolverBase::clear_assembly() {
    
    _assembly = NULL;
}




libMesh::NumericVector<Real>&
MAST::ComplexSolverBase::real_solution() {
    
    libMesh::System& sys = _assembly->system();
    
    std::string nm;
    nm = sys.name();
    nm += "real_sol";
    
    if (!sys.have_vector(nm))
        sys.add_vector(nm);

    return sys.get_vector(nm);
}



const libMesh::NumericVector<Real>&
MAST::ComplexSolverBase::real_solution() const {
    
    libMesh::System& sys = _assembly->system();

    std::string nm;
    nm = sys.name();
    nm += "real_sol";
    
    if (!sys.have_vector(nm))
        sys.add_vector(nm);
    
    return sys.get_vector(nm);
}


libMesh::NumericVector<Real>&
MAST::ComplexSolverBase::imag_solution() {
    
    libMesh::System& sys = _assembly->system();

    std::string nm;
    nm = sys.name();
    nm += "imag_sol";
    
    if (!sys.have_vector(nm))
        sys.add_vector(nm);
    
    return sys.get_vector(nm);
}


const libMesh::NumericVector<Real>&
MAST::ComplexSolverBase::imag_solution() const {
    
    libMesh::System& sys = _assembly->system();

    std::string nm;
    nm = sys.name();
    nm += "imag_sol";
    
    if (!sys.have_vector(nm))
        sys.add_vector(nm);
    
    return sys.get_vector(nm);
}



void
MAST::ComplexSolverBase::solve() {
    
    
    //  The complex system of equations
    //     (J_R + i J_I) (x_R + i x_I) - (f_R + i f_I) = 0
    //  is rewritten as
    //     [ J_R   -J_I] {x_R}  -  {f_R}  = {0}
    //     [ J_I    J_R] {x_I}  -  {f_I}  = {0}
    //
    
    
    // continue iterations till the L2 residual of both real and imaginary
    // parts is satisfied
    
    bool
    if_cont = true,
    if_re   = true;  // keeps track of whether real or imag part is being solved
    
    Real
    res_l2  = 0.;
    
    unsigned int
    iters   = 0;
    
    libMesh::System& sys = _assembly->system();

    
    while (if_cont) {
        
        // tell the assembly object to now assemble
        if (if_re) {

            // swap the solution with the real part
            *(sys.solution) =  this->real_solution();
            _assembly->set_assemble_real_part();
            libMesh::out << "Solving Real Part: " << std::endl;
        }
        else {
            
            *(sys.solution) =  this->imag_solution();
            _assembly->set_assemble_imag_part();
            libMesh::out << "Solving Imaginary Part: " << std::endl;
        }
        sys.solution->close();
        
        // tell implicit system to solve
        sys.solve();
        
        // now copy the solution vector back
        if (if_re) {
            
            this->real_solution() = (*sys.solution);
            this->real_solution().close();
            if_re = false;
        }
        else {
            
            this->imag_solution() = (*sys.solution);
            this->imag_solution().close();
            if_re = true;
        }
        
        // check the residual
        res_l2  =  _assembly->residual_l2_norm();
        libMesh::out << "Complex Residual L2-norm: " << res_l2 << std::endl;
        
        
        // increment the iteration counter
        iters++;
        
        
        if (res_l2 <= tol ||
            iters  >= max_iters) {
            
            libMesh::out
            << "Terminating complex solver iterations!" << std::endl;
            if_cont = false;
        }
    }
    
}


