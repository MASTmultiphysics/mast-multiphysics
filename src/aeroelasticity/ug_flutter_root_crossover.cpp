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


// C++ includes
#include <iomanip>


// MAST includes
#include "aeroelasticity/ug_flutter_root_crossover.h"
#include "aeroelasticity/ug_flutter_root.h"
#include "aeroelasticity/flutter_solution_base.h"


MAST::UGFlutterRootCrossover::UGFlutterRootCrossover():
MAST::FlutterRootCrossoverBase()
{ }


void
MAST::UGFlutterRootCrossover::print(std::ostream &output) const {
    
    const MAST::UGFlutterRoot
    &lower = dynamic_cast<const MAST::UGFlutterRoot&>(crossover_solutions.first->get_root(root_num)),
    &upper = dynamic_cast<const MAST::UGFlutterRoot&>(crossover_solutions.second->get_root(root_num));
    
    output
    << " Lower Root: " << std::endl
    << "    k : " << std::setw(15) << lower.kr << std::endl
    << "    g : " << std::setw(15) << lower.g << std::endl
    << "    V : " << std::setw(15) << lower.V << std::endl
    << "omega : " << std::setw(15) << lower.omega << std::endl
    << " Upper Root: " << std::endl
    << "    k : " << std::setw(15) << upper.kr << std::endl
    << "    g : " << std::setw(15) << upper.g << std::endl
    << "    V : " << std::setw(15) << upper.V << std::endl
    << "omega : " << std::setw(15) << upper.omega << std::endl;
    
    if (root) {
        const MAST::UGFlutterRoot&
        r = dynamic_cast<const MAST::UGFlutterRoot&>(*root);
        
        output
        << "Critical Root: " << std::endl
        << "    k : " << std::setw(15) << r.kr << std::endl
        << "    g : " << std::setw(15) << r.g << std::endl
        << "    V : " << std::setw(15) << r.V << std::endl
        << "omega : " << std::setw(15) << r.omega << std::endl;
    }
    else
        output << "Critical root not yet calculated." << std::endl;
}

