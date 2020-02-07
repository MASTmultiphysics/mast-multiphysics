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

//#define CATCH_CONFIG_MAIN  // This tells Catch to provide main() - only do this in one cpp file
#define CATCH_CONFIG_RUNNER
#include "catch.hpp"

#include "libmesh/libmesh.h"

libMesh::LibMeshInit *p_global_init;

int main(int argc, char* argv[])
{
    p_global_init = new libMesh::LibMeshInit(argc, (const char **)argv);
    
    int result = Catch::Session().run(argc, argv);
    
    delete p_global_init;
    
    return result;
}
