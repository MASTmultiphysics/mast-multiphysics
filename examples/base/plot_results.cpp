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


// MAST includes
#include "examples/base/plot_results.h"
#include "base/complex_mesh_field_function.h"
#include "solver/complex_solver_base.h"

// libMesh includes
#include "libmesh/numeric_vector.h"
#include "libmesh/exodusII_io.h"


void
MAST::plot_structural_flutter_solution(const std::string& nm,
                                       libMesh::System& sys,
                                       const ComplexVectorX& eig_vec,
                                       const std::vector<libMesh::NumericVector<Real>*>& basis) {
    
    // now write the flutter mode to an output file.
    // Flutter mode Y = sum_i (X_i * (xi_re + xi_im)_i)
    // using the right eigenvector of the system.
    // where i is the structural mode
    //
    // The time domain simulation assumes the temporal solution to be
    // X(t) = (Y_re + i Y_im) exp(p t)
    //      = (Y_re + i Y_im) exp(p_re t) * (cos(p_im t) + i sin(p_im t))
    //      = exp(p_re t) (Z_re + i Z_im ),
    // where Z_re = Y_re cos(p_im t) - Y_im sin(p_im t), and
    //       Z_im = Y_re sin(p_im t) + Y_im cos(p_im t).
    //
    // We write the simulation of the mode over a period of oscillation
    //
    
    
    // first calculate the real and imaginary vectors
    std::auto_ptr<libMesh::NumericVector<Real> >
    sol_copy(sys.solution->clone().release()),
    re(sys.solution->zero_clone().release()),
    im(sys.solution->zero_clone().release());
    
    
    // first the real part
    sys.solution->zero();
    for (unsigned int i=0; i<basis.size(); i++) {
        re->add(eig_vec(i).real(), *basis[i]);
        im->add(eig_vec(i).imag(), *basis[i]);
    }
    re->close();
    im->close();
    
    // now open the output processor for writing
    libMesh::ExodusII_IO flutter_mode_output(sys.get_mesh());
    
    // use N steps in a time-period
    Real
    t_sys = sys.time,
    pi    = acos(-1.);
    unsigned int
    N_divs = 100;
    
    
    for (unsigned int i=0; i<=N_divs; i++) {
        sys.time   =  2.*pi*(i*1.)/(N_divs*1.);
        
        sys.solution->zero();
        sys.solution->add( cos(sys.time), *re);
        sys.solution->add(-sin(sys.time), *im);
        sys.solution->close();
        flutter_mode_output.write_timestep(nm,
                                           sys.get_equation_systems(),
                                           i+1,
                                           sys.time);
    }
    
    // reset the system time and solution
    sys.time = t_sys;
    *sys.solution = *sol_copy;
    sys.solution->close();
}



void
MAST::plot_fluid_flutter_solution(const std::string&                      nm,
                                  libMesh::System&            structural_sys,
                                  libMesh::System&                 fluid_sys,
                                  MAST::ComplexMeshFieldFunction&      displ,
                                  MAST::ComplexSolverBase&            solver,
                                  const ComplexVectorX&               eig_vec,
                                  const std::vector<libMesh::NumericVector<Real>*>& structural_basis) {
    
    // now write the flutter mode to an output file.
    // Flutter mode Y = sum_i (X_i * (xi_re + xi_im)_i)
    // using the right eigenvector of the system.
    // where i is the structural mode
    //
    // The time domain simulation assumes the temporal solution to be
    // X(t) = (Y_re + i Y_im) exp(p t)
    //      = (Y_re + i Y_im) exp(p_re t) * (cos(p_im t) + i sin(p_im t))
    //      = exp(p_re t) (Z_re + i Z_im ),
    // where Z_re = Y_re cos(p_im t) - Y_im sin(p_im t), and
    //       Z_im = Y_re sin(p_im t) + Y_im cos(p_im t).
    //
    // We write the simulation of the mode over a period of oscillation
    //
    
    
    // next write the fluid mode
    {
        // get the size of basis from the structural node
        unsigned int
        n_basis = (unsigned int) structural_basis.size();
        
        // first calculate the fluid basis
        std::vector<libMesh::NumericVector<Real>*>
        fluid_basis_re(n_basis),
        fluid_basis_im(n_basis);
        
        std::auto_ptr<libMesh::NumericVector<Real> >
        zero(structural_sys.solution->zero_clone().release());
        
        for (unsigned int i=0; i<n_basis; i++) {
            
            // solve the fluid system for the given structural mode and
            // the frequency of the flutter root
            displ.clear();
            displ.init         (*structural_basis[i], *zero);
            
            solver.solve_block_matrix();
            
            fluid_basis_re[i] = solver.real_solution().clone().release();
            fluid_basis_im[i] = solver.imag_solution().clone().release();
        }
        
        
        // first calculate the real and imaginary vectors
        std::auto_ptr<libMesh::NumericVector<Real> >
        re(fluid_sys.solution->zero_clone().release()),
        im(fluid_sys.solution->zero_clone().release());
        
        
        // first the real part
        fluid_sys.solution->zero();
        for (unsigned int i=0; i<n_basis; i++) {
            
            re->add( eig_vec(i).real(), *fluid_basis_re[i]);
            re->add(-eig_vec(i).imag(), *fluid_basis_im[i]);
            
            im->add(eig_vec(i).imag(), *fluid_basis_re[i]);
            im->add(eig_vec(i).real(), *fluid_basis_im[i]);
        }
        re->close();
        im->close();
        
        // now open the output processor for writing
        libMesh::ExodusII_IO flutter_mode_output(fluid_sys.get_mesh());
        
        // use N steps in a time-period
        Real
        t_sys = fluid_sys.time,
        pi    = acos(-1.);
        unsigned int
        N_divs = 100;
        
        
        for (unsigned int i=0; i<=N_divs; i++) {
            fluid_sys.time   =  2.*pi*(i*1.)/(N_divs*1.);
            
            fluid_sys.solution->zero();
            fluid_sys.solution->add( cos(fluid_sys.time), *re);
            fluid_sys.solution->add(-sin(fluid_sys.time), *im);
            fluid_sys.solution->close();
            flutter_mode_output.write_timestep("fluid_flutter_mode.exo",
                                               fluid_sys.get_equation_systems(),
                                               i+1,
                                               fluid_sys.time);
        }
        
        // reset the system time
        fluid_sys.time = t_sys;
        
        // now delete the solution
        for (unsigned int i=0; i<n_basis; i++) {
            delete fluid_basis_re[i];
            delete fluid_basis_im[i];
        }
    }
}


