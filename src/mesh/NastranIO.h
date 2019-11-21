#ifndef NASTRANIO_H
#define NASTRANIO_H


// C++ Includes
#include <iostream>
#include <map>
#include <list>

// libMesh Includes
#include "libmesh/libmesh_common.h"
#include "libmesh/mesh_input.h"
#include "libmesh/elem.h"

#include "Python.h"
#include "mesh/NastranIO.h"
#include "mesh/pynastranIO.h"

void printElementMap(std::map<std::string, std::vector<std::vector<int>>> elementMap);

void printNodeCoords(std::vector<std::vector<float>> nodes);


/** 
 * The NastranIO class is a preliminary implementation for reading NASTRAN
 * mesh files using pyNastran with BDF files as input.
 */
class NastranIO : public libMesh::MeshInput<libMesh::MeshBase>
{
public:
    explicit
    NastranIO(libMesh::MeshBase& mesh, const bool pythonPreinitialized=false);
    
    virtual ~NastranIO ();
    
    virtual void read (const std::string & name) override;
    virtual void read (BDFModel* model);
    
    std::map<uint64_t, libMesh::Node*> getNastran2libMeshNodeMap();
    std::map<const libMesh::Node*, uint64_t> getlibMesh2NastranNodeMap();
    
    std::map<uint64_t, libMesh::Elem*> getNastran2libMeshElemMap();
    std::map<libMesh::Elem*, uint64_t> getlibMesh2NastranElemMap();
    
    std::map<std::pair<int,int>, int> getPIDElemtype2SubdomainIDMap();
    std::map<int, std::set<int>> getPID2subdomainIDsMap();
    
private:
    
    const bool  _pythonPreinitialized = false;
    bool        pythonInitialized =     false;
    
    std::map<uint64_t, libMesh::Node*> nastran2libMeshNodeMap;
    std::map<const libMesh::Node*, uint64_t> libMesh2NastranNodeMap;
    std::map<uint64_t, libMesh::Elem*> nastran2libMeshElemMap;
    std::map<libMesh::Elem*, uint64_t> libMesh2NastranElemMap;
    std::map<std::pair<int, int>, int> pid_elemType2subdomainMap = {};
    
    void read_nodes(BDFModel* model, libMesh::MeshBase& the_mesh);
    void read_elements(BDFModel* model, libMesh::MeshBase& the_mesh);
    
    // Map from NASTRAN Elements to Equivalent libMesh Elements
    // TODO: Not yet complete, need to add all NASTRAN elements which we need support for.
    std::map<std::string, libMesh::ElemType> nastranToLibMeshElem = {
        
        // 0D Elements (i.e. Ground Springs)
        {"CELAS1_1", libMesh::NODEELEM},    {"CELAS2_1",  libMesh::NODEELEM},
        {"CELAS3_1", libMesh::NODEELEM},    {"CELAS4_1",  libMesh::NODEELEM},
        {"CBUSH_1",  libMesh::NODEELEM},    {"CBUSH1D_1", libMesh::NODEELEM},
        
        // 1D Elements
        {"CELAS1_2", libMesh::EDGE2},       {"CELAS2_2",  libMesh::EDGE2},
        {"CELAS3_2", libMesh::EDGE2},       {"CELAS4_2",  libMesh::EDGE2},
        {"CBUSH_2",  libMesh::EDGE2},       {"CBUSH1D_2", libMesh::EDGE2},
        {"CBUSH2D_2",libMesh::EDGE2},       
        {"CROD_2",   libMesh::EDGE2},       {"CBAR_2",    libMesh::EDGE2},      
        {"CBEAM_2",  libMesh::EDGE2},       {"CBEAM3_3",  libMesh::EDGE3},
        
        // 2D Elements
        {"CTRIA3_3", libMesh::TRI3},        {"CTRIA6_6",  libMesh::TRI6},
        {"CTRIAR_3", libMesh::TRI3},
        {"CQUAD4_4", libMesh::QUAD4},       {"CQUAD8_8",  libMesh::QUAD8},
        {"CQUADR_4", libMesh::QUADSHELL4},  {"CQUAD_4",   libMesh::QUAD4},
        {"CQUAD_8",  libMesh::QUAD8},       {"CQUAD_9",   libMesh::QUAD9},
        {"CQUADX_4", libMesh::QUAD4},       {"CQUADX_8",  libMesh::QUAD8},
        {"CQUADX_9", libMesh::QUAD9},
        
        // 3D Elements
        {"CTETRA_4", libMesh::TET4},        {"CTETRA_10", libMesh::TET10},
        {"CPENTA_6", libMesh::PRISM6},      {"CPENTA_15", libMesh::PRISM15},
        {"CPYRAM_5", libMesh::PYRAMID5},    {"CPYRAM_13", libMesh::PYRAMID13},
        {"CHEXA_8",  libMesh::HEX8},        {"CHEXA_20",  libMesh::HEX20}
    };
    
    void initializePython();
    void finalizePython();
};


#endif // NASTRANIO_H
