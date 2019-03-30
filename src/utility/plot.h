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

#ifndef __mast_plot_h__
#define __mast_plot_h__

// MAST includes
#include "base/mast_data_types.h"
#include "base/mast_config.h"

// libMesh includes
#include "libmesh/elem.h"
#include "libmesh/node.h"

#if MAST_ENABLE_GNUPLOT == 1

#include "gnuplot-iostream.h"

namespace MAST {
    
    void plot_elem(Gnuplot& gp, const libMesh::Elem& elem);
    void plot_node(Gnuplot& gp, const libMesh::Node& node);
}

#endif // MAST_ENABLE_GNUPLOT

#endif // __mast_plot_h__
