/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2020  Manav Bhatia and MAST authors
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

#ifndef __mast_topology_level_set_nucleation__
#define __mast_topology_level_set_nucleation__

// MAST includes
#include "base/mast_data_types.h"


namespace MAST {

namespace Examples {

class LevelSetNucleationFunction2D:
public MAST::FieldFunction<RealVectorX> {
    
public:
    LevelSetNucleationFunction2D(Real x0,
                                 Real y0,
                                 Real l1,
                                 Real l2,
                                 Real nx_mesh,
                                 Real ny_mesh,
                                 Real nx_holes,
                                 Real ny_holes):
    MAST::FieldFunction<RealVectorX>("Phi"),
    _x0  (x0),
    _y0  (y0),
    _l1  (l1),
    _l2  (l2),
    _nx_mesh  (nx_mesh),
    _ny_mesh  (ny_mesh),
    _nx_holes (nx_holes),
    _ny_holes (ny_holes),
    _pi  (acos(-1.)) {
        
        Real
        dx = _l1/(1.*_nx_holes);
        
        for (unsigned int i=0; i<_nx_holes; i++)
            _x_axis_hole_locations.insert(_x0+(i+.5)*dx);
        
        //
        // now, along the y-axis
        //
        dx = _l2/(1.*_ny_holes);
        for (unsigned int i=0; i<_ny_holes; i++)
            _y_axis_hole_locations.insert(_y0+(i+0.5)*dx);
    }
    
    virtual ~LevelSetNucleationFunction2D() {}
    
    virtual void operator()(const libMesh::Point& p,
                            const Real t,
                            RealVectorX& v) const {
        
        libmesh_assert_less_equal(t, 1);
        libmesh_assert_equal_to(v.size(), 1);
        
        //
        // the libMesh solution projection routine for Lagrange elements
        // will query the function value at the nodes. So, we figure
        // out which nodes should have zero values set to them.
        // if there is one hole in any direction, it will be in the
        // center of the domain. If there are more than 1, then two of
        // the holes will be on the boundary and others will fill the
        // interior evenly.
        //
        const Real
        dx_mesh = _l1/(1.*_nx_holes),
        dy_mesh = _l2/(1.*_ny_holes);
        
        std::set<Real>::const_iterator
        x_it_low = _x_axis_hole_locations.lower_bound(p(0)-dx_mesh),
        y_it_low = _y_axis_hole_locations.lower_bound(p(1)-dy_mesh);
        
        unsigned int
        n = 0;
        //
        // see if the x-location needs a hole
        //
        for ( ; x_it_low != _x_axis_hole_locations.end(); x_it_low++) {
            if (std::fabs(*x_it_low - p(0)) <= dx_mesh*0.25) {
                n++;
                break;
            }
        }
        
        //
        // now check the y-location
        //
        for ( ; y_it_low != _y_axis_hole_locations.end(); y_it_low++) {
            if (std::fabs(*y_it_low - p(1)) <= dy_mesh*0.25) {
                n++;
                break;
            }
        }
        
        if (n == 2)
            v(0) = -1.e0;
        else
            v(0) = 1.e0;
    }
    
    
protected:
    Real
    _x0,
    _y0,
    _l1,
    _l2,
    _nx_mesh,
    _ny_mesh,
    _nx_holes,
    _ny_holes,
    _pi;
    std::set<Real> _x_axis_hole_locations;
    std::set<Real> _y_axis_hole_locations;
};


class LevelSetNucleationFunction3D:
public MAST::FieldFunction<RealVectorX> {
    
public:
    LevelSetNucleationFunction3D(Real x0,
                                 Real y0,
                                 Real z0,
                                 Real l1,
                                 Real l2,
                                 Real l3,
                                 Real nx_mesh,
                                 Real ny_mesh,
                                 Real nz_mesh,
                                 Real nx_holes,
                                 Real ny_holes,
                                 Real nz_holes):
    MAST::FieldFunction<RealVectorX>("Phi"),
    _x0  (x0),
    _y0  (y0),
    _z0  (z0),
    _l1  (l1),
    _l2  (l2),
    _l3  (l3),
    _nx_mesh  (nx_mesh),
    _ny_mesh  (ny_mesh),
    _nz_mesh  (nz_mesh),
    _nx_holes (nx_holes),
    _ny_holes (ny_holes),
    _nz_holes (nz_holes),
    _pi  (acos(-1.)) {
        
        Real
        dx = _l1/(1.*_nx_holes);
        
        for (unsigned int i=0; i<_nx_holes; i++)
            _x_axis_hole_locations.insert(_x0+(i+.5)*dx);
        
        //
        // now, along the y-axis
        //
        dx = _l2/(1.*_ny_holes);
        for (unsigned int i=0; i<_ny_holes; i++)
            _y_axis_hole_locations.insert(_y0+(i+0.5)*dx);

        //
        // now, along the z-axis
        //
        dx = _l3/(1.*_nz_holes);
        for (unsigned int i=0; i<_nz_holes; i++)
            _z_axis_hole_locations.insert(_z0+(i+0.5)*dx);
    }
    
    virtual ~LevelSetNucleationFunction3D() {}
    
    virtual void operator()(const libMesh::Point& p,
                            const Real t,
                            RealVectorX& v) const {
        
        libmesh_assert_less_equal(t, 1);
        libmesh_assert_equal_to(v.size(), 1);
        
        //
        // the libMesh solution projection routine for Lagrange elements
        // will query the function value at the nodes. So, we figure
        // out which nodes should have zero values set to them.
        // if there is one hole in any direction, it will be in the
        // center of the domain. If there are more than 1, then two of
        // the holes will be on the boundary and others will fill the
        // interior evenly.
        //
        const Real
        dx_mesh = _l1/(1.*_nx_holes),
        dy_mesh = _l2/(1.*_ny_holes),
        dz_mesh = _l3/(1.*_nz_holes);

        std::set<Real>::const_iterator
        x_it_low = _x_axis_hole_locations.lower_bound(p(0)-dx_mesh),
        y_it_low = _y_axis_hole_locations.lower_bound(p(1)-dy_mesh),
        z_it_low = _z_axis_hole_locations.lower_bound(p(2)-dz_mesh);

        unsigned int
        n = 0;
        //
        // see if the x-location needs a hole
        //
        for ( ; x_it_low != _x_axis_hole_locations.end(); x_it_low++) {
            if (std::fabs(*x_it_low - p(0)) <= dx_mesh*0.25) {
                n++;
                break;
            }
        }
        
        //
        // now check the y-location
        //
        for ( ; y_it_low != _y_axis_hole_locations.end(); y_it_low++) {
            if (std::fabs(*y_it_low - p(1)) <= dy_mesh*0.25) {
                n++;
                break;
            }
        }

        //
        // now check the z-location
        //
        for ( ; z_it_low != _z_axis_hole_locations.end(); z_it_low++) {
            if (std::fabs(*z_it_low - p(2)) <= dz_mesh*0.25) {
                n++;
                break;
            }
        }

        if (n == 3)
            v(0) = -1.e0;
        else
            v(0) = 1.e0;
    }
    
    
protected:
    Real
    _x0,
    _y0,
    _z0,
    _l1,
    _l2,
    _l3,
    _nx_mesh,
    _ny_mesh,
    _nz_mesh,
    _nx_holes,
    _ny_holes,
    _nz_holes,
    _pi;
    std::set<Real> _x_axis_hole_locations;
    std::set<Real> _y_axis_hole_locations;
    std::set<Real> _z_axis_hole_locations;
};

}
}


#endif // __mast_topology_level_set_nucleation__
