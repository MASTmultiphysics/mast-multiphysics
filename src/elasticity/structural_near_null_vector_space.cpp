/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2019  Manav Bhatia
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
#include "elasticity/structural_near_null_vector_space.h"
#include "base/nonlinear_system.h"

// libMesh includes
#include "libmesh/function_base.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/exodusII_io.h"

namespace MAST {
    
    class StructuralModes:
    public libMesh::FunctionBase<Real> {
    public:
        StructuralModes(const unsigned int mode_num):
        libMesh::FunctionBase<Real>(),
        _mode(mode_num) {
        
            // make sure the mode number makes sense
            libmesh_assert_less_equal(mode_num, 5);
        }
        
        
      virtual std::unique_ptr<libMesh::FunctionBase<Real> > clone () const {
            
	libMesh::FunctionBase<Real> *rval = new MAST::StructuralModes(_mode);
	return std::unique_ptr<libMesh::FunctionBase<Real> >(rval);
        }
        
        // this should not get called
        virtual Real operator()
        (const libMesh::Point& p, const Real time) {libmesh_assert(false);}
        
        virtual void
        operator() (const libMesh::Point& p,
                    const Real time,
                    libMesh::DenseVector<Real>& output) {
            
            output.resize(6);
            
            switch (_mode) {
                case 0:
                case 1:
                case 2: {
                    output(_mode) = 1.;
                }
                    break;
                    
                case 3: {
                    
                    Real th = 1.e-5;
                    // rigid-body rotation
                    output(_mode) = th;
                    // also set deformation assuming rotation about (x,y,z)=0
                    // theta-x is only going to lead to changes in y,z
                    output(1) =  (cos(th) * p(1) - sin(th) * p(2)) - p(1);
                    output(2) =  (sin(th) * p(1) + cos(th) * p(2)) - p(2);
                }
                    break;

                case 4: {
                    
                    Real th = 1.e-5;
                    // rigid-body rotation
                    output(_mode) = th;
                    // also set deformation assuming rotation about (x,y,z)=0
                    // theta-y is only going to lead to changes in x,z
                    output(0) = ( cos(th) * p(0) + sin(th) * p(2)) - p(0);
                    output(2) = (-sin(th) * p(0) + cos(th) * p(2)) - p(2);
                }
                    break;

                case 5: {
                    
                    Real th = 1.e-5;
                    // rigid-body rotation
                    output(_mode) = th;
                    // also set deformation assuming rotation about (x,y,z)=0
                    // theta-z is only going to lead to changes in x,y
                    output(0) =  (cos(th) * p(0) - sin(th) * p(1)) - p(0);
                    output(1) =  (sin(th) * p(0) + cos(th) * p(1)) - p(1);
                }
                    break;

                default:
                    libmesh_error(); // should not get here
            }
        }
    protected:
        
        /*!
         *  mode number: 
         *  0 = rigid-body translation in x
         *  1 = rigid-body translation in y
         *  2 = rigid-body translation in z
         *
         *  3 = rigid-body rotation about x
         *  4 = rigid-body rotation about y
         *  5 = rigid-body rotation about z
         */
        const unsigned int _mode;
    };
}



MAST::StructuralNearNullVectorSpace::StructuralNearNullVectorSpace():
libMesh::NonlinearImplicitSystem::ComputeVectorSubspace() {

    _nm.resize(6);
    _nm[0]   =  "x_translation_mode";
    _nm[1]   =  "y_translation_mode";
    _nm[2]   =  "z_translation_mode";
    _nm[3]   =  "tx_rotation_mode";
    _nm[4]   =  "ty_rotation_mode";
    _nm[5]   =  "tz_rotation_mode";
    
}


void
MAST::StructuralNearNullVectorSpace::
operator()(std::vector<libMesh::NumericVector<Real>*>&sp,
           libMesh::NonlinearImplicitSystem& s) {
    

    MAST::NonlinearSystem& sys = dynamic_cast<MAST::NonlinearSystem&>(s);
    
    // make sure that the vector coming in is of zero-size
    libmesh_assert_equal_to(sp.size(), 0);
    
    // now initialize six vectors and add them to the system and sp
    sp.resize(6);
    
    // the names of these six vectors, for storage
    std::vector<std::string> nm(6);
    nm[0]   =  "x_translation_mode";
    nm[1]   =  "y_translation_mode";
    nm[2]   =  "z_translation_mode";
    nm[3]   =  "tx_rotation_mode";
    nm[4]   =  "ty_rotation_mode";
    nm[5]   =  "tz_rotation_mode";
    
    
    // first, the three translation modes
    for (unsigned int i=0; i<nm.size(); i++) {
        
        sp[i] = &s.add_vector(nm[i]);
        sp[i]->zero();
        
        MAST::StructuralModes sol(i);
        sys.project_vector_without_dirichlet(*sp[i], sol);
        
        //libMesh::ExodusII_IO(s.get_mesh()).write_equation_systems("output.exo",
        //                                                          s.get_equation_systems());
        
    }
}


