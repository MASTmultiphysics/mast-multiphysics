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

#ifndef __mast__mesh_initializer_h__
#define __mast__mesh_initializer_h__


// C++ includes
#include <vector>
#include <map>

// MAST includes
#include "base/mast_data_types.h"

// libMesh includes
#include "libmesh/libmesh_config.h"
#include "libmesh/unstructured_mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/node.h"

namespace MAST {
    
    
    
    /*!
     *   This class distributes mesh points over a mesh. Input requires definition
     *   of divisions along each mesh dimension.
     */
    
    class MeshInitializer {
        
    public:
        
        MeshInitializer():
        _mesh(nullptr)
        { }
        
        virtual ~MeshInitializer()
        { }
        
        /*!
         *   This class defies the divisions along a computational coordinate
         *   of the mesh.
         */
        class CoordinateDivisions{
        protected:
            
            /*!
             *   number of divisions along this coordinate
             */
            unsigned int _n_divs;
            
            /*!
             *   location at each division along the coordinate.
             */
            std::vector<Real> _div_locations;
            
            /*!
             *   number of elements within each division
             */
            std::vector<unsigned int> _n_subdivs_in_div;
            
            /*!
             *   relative numbers specifying the expected mesh size at each
             *   division point
             */
            std::vector<Real> _relative_mesh_size_at_div;
            
            /*!
             *   map of computational and nodal location along the coordinate
             */
            std::vector<Real> _points;
            
        public:
            
            CoordinateDivisions():
            _n_divs(0)
            { }
            
            void init(unsigned int n, const std::vector<Real>& div_loc,
                      const std::vector<Real>& dx_relative,
                      const std::vector<unsigned int>& n_dx)
            {
                libmesh_assert ( n > 0 );
                libmesh_assert ( div_loc.size()     == n+1 );
                libmesh_assert ( n_dx.size()        == n   );
                libmesh_assert ( dx_relative.size() == n+1 );
                
                
                _n_divs = n;
                
                _div_locations = div_loc;
                _n_subdivs_in_div = n_dx;
                _relative_mesh_size_at_div = dx_relative;
                _points.clear();
                this->distributePoints();
            }
            
            
            /*!
             *   @returns number of division along this coordinate
             */
            unsigned int n_divs() const {
                libmesh_assert(_n_divs > 0);
                return _n_divs;
            }
            
            /*!
             *   @returns location of division point \p i
             */
            Real div_location(unsigned int i) const {
                libmesh_assert(_n_divs > 0);
                libmesh_assert(i <= _n_divs);
                return _div_locations[i];
            }

            
            /*!
             *   @returns relative mesh size at point \p i
             */
            Real relative_mesh_size(unsigned int i) const {
                libmesh_assert(_n_divs > 0);
                libmesh_assert(i <= _n_divs);
                return _relative_mesh_size_at_div[i];
            }

            
            /*!
             *   @returns number of elements in division \p i
             */
            Real n_elements_in_div(unsigned int i) const {
                libmesh_assert(_n_divs > 0);
                libmesh_assert(i < _n_divs);
                return _n_subdivs_in_div[i];
            }
            
            
            /*!
             *    @returns the length of the mesh along this coordinate
             */
            Real length() const {
                libmesh_assert( _n_divs > 0);
                return (*_div_locations.rbegin());
            }
            
            
            /*!
             *    @returns the length of the mesh along this coordinate
             */
            unsigned int total_elem_divs() const {
                libmesh_assert( _n_divs > 0);
                unsigned int n=0;
                std::vector<unsigned int>::const_iterator it = _n_subdivs_in_div.begin();
                for ( ; it != _n_subdivs_in_div.end(); it++)
                    n += *it;
                return n;
            }
            
            /*!
             *   @returns the location of the specified node
             */
            Real operator() (const Real eta) {
                libmesh_assert ((eta >= 0.) && (eta <= 1.));
                
                // idx provides the approximate location of the computational coordinate
                // in the calculate point distribution
                unsigned int idx = floor(eta * (_points.size()-1));
                
                if (eta == 0.)
                    return (*_points.begin());
                else if (eta == 1.)
                    return (*_points.rbegin());
                else
                {
                    Real idx_eta0 = Real(idx)/Real(_points.size()-1),
                    idx_eta1 = Real(idx+1)/Real(_points.size()-1);
                    return _points[idx] + (eta-idx_eta0)/(idx_eta1-idx_eta0)*
                    (_points[idx+1]-_points[idx]);
                }
            }
            
            /*!
             *   creates the nodal locations using the provided input
             */
            void distributePoints()
            {
                // make sure sizes are correct
                libmesh_assert ( _n_divs > 0 );
                
                // calculate total number of points
                unsigned int n_total_points = this->total_elem_divs() + 1;
                
                // resize the points vector and set the first and last points of each division
                _points.resize(n_total_points);
                
                unsigned int n=1;
                _points[0] = _div_locations[0];
                for (unsigned int i=0; i<_n_divs; i++)
                {
                    n += _n_subdivs_in_div[i];
                    _points[n-1] = _div_locations[i+1];
                }
                
                n=1;
                Real dx=0.0, growth_factor = 0.0;
                // now calculate the base mesh size, and calculate the nodal points
                for (unsigned int i=0; i<_n_divs; i++)
                {
                    growth_factor =
                    pow(_relative_mesh_size_at_div[i+1]/_relative_mesh_size_at_div[i],
                        1.0/(_n_subdivs_in_div[i]-1.0));
                    if (fabs(growth_factor-1.0)>1.0e-10)
                        dx = (_div_locations[i+1]-_div_locations[i]) *
                        (1.0-growth_factor)/(1.0-pow(growth_factor, _n_subdivs_in_div[i]));
                    else
                    {
                        growth_factor = 1.0;
                        dx = (_div_locations[i+1]-_div_locations[i]) / _n_subdivs_in_div[i];
                    }
                    
                    for (unsigned int n_pt=1; n_pt<_n_subdivs_in_div[i]; n_pt++)
                    {
                        _points[n+n_pt-1] = _points[n+n_pt-2] + dx;
                        dx *= growth_factor;
                    }
                    n += _n_subdivs_in_div[i];
                }
            }
            
        };
        
        /*!
         *   initializes the object with the division for each dimension. Sets-up the
         *   mesh for the provided information, and then uses the provided funciton
         *   move the mesh points.
         */
        void init (const std::vector<MeshInitializer::CoordinateDivisions*>& divs,
                   libMesh::UnstructuredMesh& mesh, libMesh::ElemType t);
        
        
    protected:
        
        /*!
         *    function modifies the mesh and sets boudnary conditions. This needs to
         *    be implemented for each inherited class
         */
        virtual void process_mesh() { }
        
        /*!
         *    mesh associated with this initialization object
         */
        libMesh::UnstructuredMesh* _mesh;
    };
    
    
}



#endif  //__MAST_point_distribution_h__
