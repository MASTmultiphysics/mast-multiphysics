//
//  interface.cpp
//  PythonC
//
//  Created by Manav Bhatia on 6/21/17.
//  Copyright © 2017 Manav Bhatia. All rights reserved.
//

#include "mast_interface.h"

void add_node(int i, double x, double y, double z, void* d) {
    
    std::cout
    << "Node: " << i
    << " x: " << x
    << " y: " << y
    << " z: " << z
    << " d*: " << d
    << std::endl;
}

void add_nothing(void *d) {
    
    std::cout
    << "Nothing"
    << " d*: " << d
    << std::endl;
}

