
// MAST includes
#include "base/system_initialization.h"


MAST::SystemInitialization::SystemInitialization (libMesh::System& sys,
                                                  const std::string& prefix):
_system(sys),
_prefix(prefix)
{ }



MAST::SystemInitialization::~SystemInitialization()
{ }

