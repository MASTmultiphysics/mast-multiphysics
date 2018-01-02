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

// C++ includes
#include <iomanip>


// MAST includes
#include "aeroelasticity/time_domain_flutter_solution.h"
#include "aeroelasticity/time_domain_flutter_root.h"
#include "numerics/lapack_dggev_interface.h"


MAST::TimeDomainFlutterSolution::TimeDomainFlutterSolution():
MAST::FlutterSolutionBase()
{ }



MAST::TimeDomainFlutterSolution::~TimeDomainFlutterSolution()
{ }




void
MAST::TimeDomainFlutterSolution::init (const MAST::TimeDomainFlutterSolver& solver,
                                       const Real v_ref,
                                       const MAST::LAPACK_DGGEV& eig_sol) {
    
    // make sure that it hasn't already been initialized
    libmesh_assert(!_roots.size());
    
    _ref_val           = v_ref;
    
    // iterate over the roots and initialize the vector
    _Amat              = eig_sol.A();
    _Bmat              = eig_sol.B();
    const ComplexMatrixX
    &VR                = eig_sol.right_eigenvectors(),
    &VL                = eig_sol.left_eigenvectors();
    const ComplexVectorX
    &num               = eig_sol.alphas();
    
    const RealVectorX
    &den               = eig_sol.betas();
    
    // use only half the modes that have a positive frequency.
    unsigned int nvals = (int)_Bmat.rows();
    
    _roots.resize(nvals);
    for (unsigned int i=0; i<nvals; i++) {
        
        MAST::TimeDomainFlutterRoot* root = new MAST::TimeDomainFlutterRoot;
        root->init(v_ref,
                   num(i),
                   den(i),
                   _Bmat,
                   VR.col(i),
                   VL.col(i));
        
        _roots[i] = root;
    }
}



unsigned int
MAST::TimeDomainFlutterSolution::
n_unstable_roots_in_upper_complex_half (Real tol) const {
    
    unsigned int n=0;
    
    // iterate over all the roots and find the ones with g >= 0
    std::vector<MAST::FlutterRootBase*>::const_iterator
    it   = _roots.begin(),
    end  = _roots.end();
    
    for ( ; it != end; it++) {
        
        // only look at upper half of the complex domain
        if ((**it).root.imag() >= 0. &&
            (**it).root.real() >= tol)
            n++;
    }
    
    return n;
}





MAST::FlutterRootBase*
MAST::TimeDomainFlutterSolution::get_critical_root(Real tol) {
    
    // If there is an unstable root, then find one with the lowest velocity
    // and send it back.
    //
    // Otherwise, find the root that is closest to the g=0 condition
    
    MAST::FlutterRootBase
    *unstable  = nullptr,
    *max_g     = nullptr;
    
    // iterate over all the roots and find the ones with g >= 0
    std::vector<MAST::FlutterRootBase*>::const_iterator
    it   = _roots.begin(),
    end  = _roots.end();
    
    for ( ; it != end; it++) {
        // only look at the upper half of the complex domain.
        if ((**it).root.imag() >= 0.) {
            
            // look for the lowest g root
            if (!max_g)
                max_g = *it;
            else if ((**it).root.real() >= max_g->root.real())
                max_g = *it;
            
            // now look for the unstable root with lowest velocity
            // only roots evaluated at V > 0 are considered, since in the absensce
            // of aerodynamics and structural damping, the roots would all
            // have zero damping.
            
            if ((**it).root.real() >= tol &&
                (**it).V            > 0.) {
                
                if (!unstable)
                    unstable = *it;
                else if ( (**it).V <= unstable->V )
                    unstable = *it;
            }
        }
        
    }
    
    if (unstable)
        return unstable;
    else
        return max_g;
}






void
MAST::TimeDomainFlutterSolution::sort(const MAST::FlutterSolutionBase& sol) {
    
    const unsigned int nvals = this->n_roots();
    libmesh_assert_equal_to(nvals, sol.n_roots());
    
    // two roots with highest modal_participation dot product are sorted
    // in the same serial order
    for (unsigned int i=0; i<nvals-1; i++)
    {
        const MAST::FlutterRootBase& r = sol.get_root(i);
        Real max_val = 0.;
        Complex val = 0.;
        unsigned int max_val_root = nvals-1;
        for (unsigned int j=i; j<nvals; j++) {
            
            // use a combination of both the eigenvectors from both
            // roots
            val = .5*(r.eig_vec_left.dot(_Bmat*_roots[j]->eig_vec_right) +
                      _roots[j]->eig_vec_left.dot(_Bmat*r.eig_vec_right));
            //_roots[j]->modal_participation.dot(r.modal_participation);
            // scale by the eigenvalue separation with the assumption that
            // the roots will be closer to each other than any other
            // root at two consecutive eigenvalues. In other words,
            // we are penalizing the dot product with the eigenvalue
            // distance
            val /= abs(r.root-_roots[j]->root);
            if (abs(val) > max_val) {
                max_val = abs(val);
                max_val_root = j;
            }
        }
        
        // now we should have the one with highest dot product
        if (i != max_val_root)
            std::swap(_roots[i], _roots[max_val_root]);
    }
}






void
MAST::TimeDomainFlutterSolution::print(std::ostream &output) {
    
    // make sure that some roots exist for this solution
    libmesh_assert(this->n_roots() > 0);
    
    const unsigned int nvals = this->n_roots();
    unsigned int n_participation_vals =
    (unsigned int)this->get_root(0).modal_participation.size();
    libmesh_assert(nvals);
    
    // first write the reference values of the root
    output << " Flutter Root " << std::endl;
    
    // now write the root
    output
    << std::setw(5) << "#"
    << std::setw(15) << "V_ref"
    << std::setw(15) << "Re"
    << std::setw(15) << "Im"
    << std::setw(3)  << " | ";
    
    // output the headers for flutter mode participation
    for (unsigned int i=0; i<n_participation_vals; i++)
        output
        << std::setw(2) << " "
        << std::setw(5) << "Mode "
        << std::setw(5) << i
        << std::setw(2) << " ";
    
    // output the headers for flutter mode
    for (unsigned int i=0; i<nvals; i++)
        output
        << std::setw(10) << "|         "
        << std::setw(5) << "Mode "
        << std::setw(5) << i
        << std::setw(10) << "         |";
    output << std::endl;
    
    for (unsigned int i=0; i<nvals; i++)
    {
        const MAST::FlutterRootBase& root = this->get_root(i);
        
        // flutter root details
        output
        << std::setw(5) << i
        << std::setw(15) << root.V
        << std::setw(15) << std::real(root.root)
        << std::setw(15) << std::imag(root.root)
        << std::setw(3)  << " | ";
        
        // now write the modal participation
        for (unsigned int j=0; j<n_participation_vals; j++)
            output
            << std::setw(12) << root.modal_participation(j)
            << std::setw(2) << " ";
        
        // now write the flutter mode
        for (unsigned int j=0; j<nvals; j++)
        {
            output
            << std::setw(2) << "| "
            << std::setw(12) << std::real(root.eig_vec_right(j))
            << std::setw(2) << " "
            << std::setw(12) << std::imag(root.eig_vec_right(j))
            << std::setw(2) << " |";
        }
        output << std::endl;
    }
}


