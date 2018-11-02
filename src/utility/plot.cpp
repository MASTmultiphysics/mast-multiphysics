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

// MAST includes
#include "utility/plot.h"

#if MAST_ENABLE_MATPLOTLIB == 1

void
MAST::plot_elem(const libMesh::Elem& elem) {
    
    // currently only for 2D.
    libmesh_assert_equal_to(elem.dim(), 2);

    unsigned int
    n_nodes = elem.n_nodes();

    std::vector<Real>
    x(n_nodes+1, 0.),
    y(n_nodes+1, 0.);

    for (unsigned int i=0; i<=n_nodes; i++) {
        x[i] = (*elem.node_ptr(i%n_nodes))(0);
        y[i] = (*elem.node_ptr(i%n_nodes))(1);
    }

    plt::plot(x, y);
    plt::pause(0.01);
}


void
MAST::plot_node(const libMesh::Node& node) {

    std::vector<Real>
    x (1, node(0)),
    y (1, node(1));
    
    plt::plot(x, y, "*");
    plt::pause(0.01);
}

#endif // MAST_ENABLE_MATPLOTLIB == 1

