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
