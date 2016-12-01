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
#include <iomanip>

// MAST includes
#include "aeroelasticity/pk_flutter_root_crossover.h"
#include "aeroelasticity/flutter_root_base.h"
#include "aeroelasticity/flutter_solution_base.h"



MAST::PKFlutterRootCrossover::PKFlutterRootCrossover():
MAST::FlutterRootCrossoverBase()
{ }


void
MAST::PKFlutterRootCrossover::print(std::ostream &output) const {
    
    const MAST::FlutterRootBase
    &lower = crossover_solutions.first->get_root(root_num),
    &upper = crossover_solutions.second->get_root(root_num);
    
    output
    << " Lower Root: " << std::endl
    << "    V : " << std::setw(15) << lower.V << std::endl
    << "   kr : " << std::setw(15) << lower.kr << std::endl
    << "   Re : " << std::setw(15) << std::real(lower.root) << std::endl
    << "   Im : " << std::setw(15) << std::imag(lower.omega) << std::endl
    << " Upper Root: " << std::endl
    << "    V : " << std::setw(15) << upper.V << std::endl
    << "   kr : " << std::setw(15) << upper.kr << std::endl
    << "   Re : " << std::setw(15) << std::real(upper.root) << std::endl
    << "   Im : " << std::setw(15) << std::imag(upper.root) << std::endl;
    
    if (root)
        output
        << "Critical Root: " << std::endl
        << "    V : " << std::setw(15) << root->V << std::endl
        << "   kr : " << std::setw(15) << root->kr << std::endl
        << "   Re : " << std::setw(15) << std::real(root->root) << std::endl
        << "   Im : " << std::setw(15) << std::imag(root->root) << std::endl;
    else
        output << "Critical root not yet calculated." << std::endl;
}
