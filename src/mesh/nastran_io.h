#ifndef __mast_nastran_io_h__
#define __mast_nastran_io_h__

// C++ includes.
#include <iostream>
#include <map>
#include <list>

// libMesh includes.
#include <libmesh/libmesh_common.h>
#include <libmesh/mesh_input.h>
#include <libmesh/elem.h>

// Python includes.
#include <Python.h>

// MAST includes.
#include "mesh/nastran_io.h"
#include "mesh/pynastran_io.h"


namespace MAST {

    void printElementMap(std::map<std::string, std::vector<std::vector<int>>> elementMap);

    void printNodeCoords(std::vector<std::vector<double>> nodes);

    /**
     * Nastran BDF mesh input.
     * The NastranIO class is a preliminary implementation for reading NASTRAN mesh information
     * using pyNastran with BDF data as input. Currently nodes, elements, & node boundary condition
     * definition are supported.
     *
     * Nodal BCs are mapped into libMesh node boundary sets based on
     * SPC IDs and nodes assigned to them in the BDF file.
     */
    class NastranIO : public libMesh::MeshInput<libMesh::MeshBase>
    {
    public:
        /**
         * Constructor.
         * @param mesh a libMesh mesh object.
         * @param python_preinit bool describing if Python has been already initialized somewhere
         *                       in the current program (by another C++/Python interface).
         */
        explicit NastranIO(libMesh::MeshBase& mesh, const bool python_preinit=false);

        /**
         * Destructor.
         */
        virtual ~NastranIO ();

        virtual void read (const std::string & name) override;
        virtual void read (BDFModel* model);

        std::map<uint64_t, libMesh::Node*> getNastran2libMeshNodeMap();
        std::map<const libMesh::Node*, uint64_t> getlibMesh2NastranNodeMap();

        std::map<uint64_t, libMesh::Elem*> getNastran2libMeshElemMap();
        std::map<libMesh::Elem*, uint64_t> getlibMesh2NastranElemMap();

        std::map<std::pair<int,int>, int> getPIDElemtype2SubdomainIDMap();
        std::map<int, std::set<int>> getPID2subdomainIDsMap();

        /**
         * Print map between Nastran property ID's (PID) to libMesh subdomain
         * ID's (SID) to libMesh::out. Note that some PID will correspond to
         * multiple SID since libMesh requires all elements in a subdomain to
         * be the same type, but Nastran allows property assignment to multiple
         * element types from the same property card (ie. one PSHELL card
         * providing properties to both CQUAD4 and CTRIA3).
         */
        void print_pid_to_subdomain_id_map();

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
        void read_node_boundaries(BDFModel* model, libMesh::MeshBase& the_mesh);

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
}

#endif // __mast_nastran_io_h__
