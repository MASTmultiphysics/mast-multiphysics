
// MAST includes
#include "base/system_initialization.h"

// libMesh includes
#include "libmesh/mesh_base.h"


MAST::SystemInitialization::SystemInitialization (libMesh::System& sys,
                                                  const std::string& prefix):
_system(sys),
_prefix(prefix) {

    // initialize the point locator for this mesh
    sys.system().get_mesh().sub_point_locator();
}



MAST::SystemInitialization::~SystemInitialization()
{ }

