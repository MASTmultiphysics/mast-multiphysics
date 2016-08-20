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

// C++ includes
#include <iostream>
#include <fstream>
#include <iomanip>
#include <limits>


// MAST includes
#include "examples/fsi/base/gaf_database.h"
#include "base/parameter.h"
#include "aeroelasticity/frequency_function.h"


MAST::GAFDatabase::GAFDatabase(const unsigned int n_modes):
MAST::FSIGeneralizedAeroForceAssembly(),
_if_evaluate(true),
_n_modes(n_modes) { }



void
MAST::GAFDatabase::set_evaluate_mode(bool f) {
    
    _if_evaluate = f;
}


void
MAST::GAFDatabase::write_gaf_file(const std::string& nm,
                                  std::vector<libMesh::NumericVector<Real>*>& modes) {
    
    libMesh::out
    << " **** Writing GAF database to : " << nm
    << "   ....  ";
    
    std::ofstream out;
    
    const unsigned int
    my_rank = modes[0]->comm().rank();
    
    if (my_rank == 0) {
        out.open(nm, std::ofstream::out);
    }
    
    const unsigned int
    n_vec_dofs = modes[0]->size(),
    n_kr       = _kr_to_gaf_map.size();
    
    
    if (my_rank == 0) {
        
        out
        << _n_modes   << std::endl
        << n_vec_dofs << std::endl
        << n_kr       << std::endl;
        
        out.precision(std::numeric_limits<Real>::max_digits10);
        
        // write the GAFs at each kr first
        {
            std::map<Real, ComplexMatrixX>::const_iterator
            it  = _kr_to_gaf_map.begin(),
            end = _kr_to_gaf_map.end();
            
            for ( ; it != end; it++) {
                
                out << it->first << std::endl;
                
                const ComplexMatrixX& mat = it->second;
                libmesh_assert(mat.rows() == _n_modes);
                libmesh_assert(mat.cols() == _n_modes);
                
                for (unsigned int i=0; i<_n_modes; i++) {
                    for (unsigned int j=0; j<_n_modes; j++)
                        out
                        << std::setw(35) << mat(i,j).real()
                        << std::setw(35) << mat(i,j).imag();
                    out << std::endl;
                }
            }
        }
        
        
        // next, write the GAF sensitivity wrt kr at each kr
        {
            std::map<Real, ComplexMatrixX>::const_iterator
            it  = _kr_to_gaf_kr_sens_map.begin(),
            end = _kr_to_gaf_kr_sens_map.end();
            
            for ( ; it != end; it++) {
                
                out << it->first << std::endl;
                
                const ComplexMatrixX& mat = it->second;
                libmesh_assert_equal_to(mat.rows(), _n_modes);
                libmesh_assert_equal_to(mat.cols(), _n_modes);
                
                for (unsigned int i=0; i<_n_modes; i++) {
                    for (unsigned int j=0; j<_n_modes; j++)
                        out
                        << std::setw(35) << mat(i,j).real()
                        << std::setw(35) << mat(i,j).imag();
                    out << std::endl;
                }
            }
        }
    }
    
    // now write the modes
    std::vector<Real> mode_vec(n_vec_dofs, 0.);
    
    for (unsigned int i=0; i<_n_modes; i++) {
        libMesh::NumericVector<Real>&
        vec = *modes[i];
        
        libmesh_assert_equal_to(vec.size(), n_vec_dofs);
        std::fill(mode_vec.begin(), mode_vec.end(), 0.);
        vec.localize(mode_vec);
        
        if (my_rank == 0) {
            
            for (unsigned int i=0; i<n_vec_dofs; i++)
                out
                << std::setw(35) << mode_vec[i] << std::endl;
            out << std::endl;
        }
    }
    
    libMesh::out
    << "   Done! " << std::endl;
    
}



void
MAST::GAFDatabase::read_gaf_file(const std::string& nm,
                                 std::vector<libMesh::NumericVector<Real>*>& modes) {
    
    libMesh::out
    << " **** Reading GAF database from : " << nm
    << "   ....  ";
    
    std::ifstream input;
    input.open(nm, std::ofstream::in);
    
    unsigned int
    n_vec_dofs = 0,
    n_kr       = 0;
    
    input >> _n_modes >> n_vec_dofs >> n_kr;
    
    Real
    kr = 0.,
    re = 0.,
    im = 0.;
    
    ComplexMatrixX
    mat (ComplexMatrixX::Zero(_n_modes, _n_modes));
    
    
    // read the GAFs at each kr first
    for (unsigned int i=0; i<n_kr; i++) {
        
        mat.setZero();
        
        input >> kr;
        
        for (unsigned int j=0; j<_n_modes; j++)
            for (unsigned int k=0; k<_n_modes; k++) {
                
                input >> re >> im;
                mat(j,k) = std::complex<Real>(re, im);
            }
        
        _kr_to_gaf_map[kr] = mat;
    }
    
    
    // next, read the GAF sensitivity wrt kr at each kr
    for (unsigned int i=0; i<n_kr; i++) {
        
        mat.setZero();
        
        input >> kr;
        
        for (unsigned int j=0; j<_n_modes; j++)
            for (unsigned int k=0; k<_n_modes; k++) {
                
                input >> re >> im;
                mat(j,k) = std::complex<Real>(re, im);
            }
        
        _kr_to_gaf_kr_sens_map[kr] = mat;
    }
    
    
    // now write the modes
    {
        std::vector<Real> mode_vec(n_vec_dofs, 0.);
        
        for (unsigned int i=0; i<_n_modes; i++) {
            libMesh::NumericVector<Real>&
            vec = *modes[i];
            
            libmesh_assert(vec.size() == n_vec_dofs);
            std::fill(mode_vec.begin(), mode_vec.end(), 0.);
            
            for (unsigned int i=0; i<n_vec_dofs; i++)
                input >> mode_vec[i];
            
            vec = mode_vec;
        }
    }
    
    libMesh::out
    << "   Done! " << std::endl;
}


ComplexMatrixX&
MAST::GAFDatabase::add_kr_mat(const Real kr,
                              const ComplexMatrixX& mat,
                              const bool if_kr_sens) {
    
    libmesh_assert_equal_to(mat.rows(), _n_modes);
    libmesh_assert_equal_to(mat.cols(), _n_modes);
    
    if (!if_kr_sens) {
        
        _kr_to_gaf_map[kr] = mat;
        return _kr_to_gaf_map[kr];
    }
    else {
        
        _kr_to_gaf_kr_sens_map[kr] = mat;
        return _kr_to_gaf_kr_sens_map[kr];
    }
}


ComplexMatrixX
MAST::GAFDatabase::get_kr_mat(const Real kr,
                              const std::map<Real, ComplexMatrixX>& data) {
    
    //
    // the following is used for calculation of the return value
    //   f(x) is defined for x for each x0 < x < x1
    //   if   x <= x0,      f(x) = f(x0)
    //   if   x0 < x < x1,  f(x) is interpolated
    //   if   x >= x1,      f(x) = f(x1)
    //
    
    ComplexMatrixX
    mat = ComplexMatrixX::Zero(_n_modes, _n_modes);
    
    std::map<Real, ComplexMatrixX>::const_iterator
    it1, it2;
    std::map<Real, ComplexMatrixX>::const_reverse_iterator
    rit  = data.rbegin();
    it1  = data.begin();
    
    // check the lower bound
    if (kr <=  it1->first) {
        mat = it1->second;
    }
    // check the upper bound
    else if (kr >=  rit->first) {
        mat = rit->second;
    }
    else {
        // if it gets here, the ordinate is in between the provided range
        it2 = data.lower_bound(kr);
        // this cannot be the first element of the map
        libmesh_assert(it2 != data.begin());
        // it2 provides the upper bound. The lower bound is provided by the
        // preceding iterator
        it1 = it2--;
        // now interpolate
        mat =  it1->second +
        (kr - it1->first)/(it2->first - it1->first) *
        (it2->second - it1->second);
    }
    
    return mat;
}


void
MAST::GAFDatabase::assemble_generalized_aerodynamic_force_matrix
(std::vector<libMesh::NumericVector<Real>*>& basis,
 ComplexMatrixX& mat,
 MAST::Parameter* p) {
    
    if (_if_evaluate) {
        
        MAST::FSIGeneralizedAeroForceAssembly::
        assemble_generalized_aerodynamic_force_matrix(basis, mat, p);
    }
    else {
        
        Real kr = 0.;
        (*_freq)(kr);
        if (!p)
            mat = this->get_kr_mat(kr, _kr_to_gaf_map);
        else
            mat = this->get_kr_mat(kr, _kr_to_gaf_kr_sens_map);
    }
}


