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

#include "base/mast_config.h"
#if MAST_ENABLE_NASTRANIO == 1

// Python includes.
#include <Python.h>

// MAST includes.
#include "mesh/nastran_io.h"
#include "mesh/pynastran_io.h"


namespace MAST {

/**
 * Nastran BDF mesh input.
 * The NastranIO class is a preliminary implementation for reading NASTRAN mesh information
 * using pyNastran with BDF data as input. Note that this class is not designed to be a complete
 * Nastran BDF reader with solution & case control, but rather a way to get basic mesh data defined
 * in BDF format into libMesh/MAST. We define basic mesh data as:
 *   - *nodes*: Nastran grids (only supports grid definition in global coordinate system)
 *   - *elements*: Nastran elements with nodal connectivity (property IDs mapped to subdomains)
 *   - *subdomains*: libMesh mesh subdomains are similar to property card assignment in Nastran
 *                   (used to connect properties to elements)
 *   - *node boundary domains*: similar to SPC ID sets in Nastran. We don't use actual BC values
 *                              assigned on SPC cards, but rather track which nodes are used in each
 *                              SPC ID. These become node boundary domains in libMesh/MAST, to which
 *                              different boundary conditions can be assigned.
 *
 * TODO: Unit tests for NastranIO class.
 */
class NastranIO : public libMesh::MeshInput<libMesh::MeshBase> {
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
        virtual ~NastranIO();

        /**
         * Read ASCII NASTRAN BDF given by filename.
         * @param filename string path to Nastran BDF formatted file.
         */
        virtual void read(const std::string & filename) override;

        /**
         * Read data directly from BDFModel object.
         * @param model pointer to BDFModel object.
         */
        virtual void read(BDFModel* model);

        /**
         * Returns mapping between Nastran BDF grid ID's and libMesh/MAST nodes.
         * @return map between Nastran grid ID's and pointers to corresponding libMesh/MAST nodes.
         */
        std::map<uint64_t, libMesh::Node*> get_nastran_to_libmesh_node_map();

        /**
         * Returns mapping between libMesh/MAST nodes and Nastran BDF grid ID's.
         * @return map between pointers to libMesh/MAST nodes and corresponding Nastran grid ID's.
         */
        std::map<const libMesh::Node*, uint64_t> get_libmesh_to_nastran_node_map();

        /**
         * Returns mapping between Nastran BDF element ID's and libMesh/MAST elements.
         * @return map between Nastran element ID's and pointers to corresponding
         *         libMesh/MAST elements.
         */
        std::map<uint64_t, libMesh::Elem*> get_nastran_to_libmesh_elem_map();

        /**
         * Returns mapping between libMesh/MAST elements and Nastran BDF element ID's.
         * @return map between pointer to libMesh/MAST elements and corresponding Nastran
         *         element ID's.
         */
        std::map<libMesh::Elem*, uint64_t> get_libmesh_to_nastran_elem_map();

        /**
         * Provides mapping between Nastran property ID & element-type pairs and corresponding
         * libMesh/MAST subdomain ID's. Note we must use a pair combining a Nastran property ID
         * along with a specific element type to get a unique map key to libMesh/MAST subdomains.
         * This is because Nastran allows one property ID to be assigned to elements of different
         * type (ie. PSHELL to both CQUAD4 & CTRIA3), but libMesh does not allow elements of
         * different type to belong to the same subdomain.
         * @return map between Nastran property ID + element type pairs and corresponding
         *         libMesh/MAST subdomain ID's.
         */
        std::map<std::pair<int,int>, int> get_nastran_pid_elemtype_to_libmesh_subdomain_map();

        /**
         * Returns mapping between Nastran property ID's and sets of libMesh/MAST subdomains. Note
         * that one Nastran property ID can map to multiple libMesh/MAST subdomains because Nastran
         * allows one property ID to be assigned to elements of different type (ie. PSHELL to both
         * CQUAD4 & CTRIA3), but libMesh does not allow elements of different type to belong to
         * the same subdomain. In this case, we use multiple subdomains to contain the different
         * element types.
         * @return map between Nastran property ID's and sets of corresponding libMesh/MAST
         *         subdomain ID's.
         */
        std::map<int, std::set<int>> get_nastran_property_to_libmesh_subdomains_map();

        /**
         * Print map between Nastran property ID's (PID) to libMesh subdomain ID's (SID) to
         * libMesh::out. Note that some PID will correspond to multiple SID since libMesh requires
         * all elements in a subdomain to be the same type, but Nastran allows property assignment
         * to multiple element types from the same property card (ie. PSHELL to both CQUAD4 and
         * CTRIA3).
         */
        void print_pid_to_subdomain_id_map();

    private:

        /// Indicates if Python was initialized outside of NastranIO class.
        const bool  python_preinitialized = false;
        /// Indicates is Python has been initialized.
        bool        python_initialized =     false;

        /// Mapping from Nastran grid IDs from BDF input to pointers to libMesh/MAST nodes.
        std::map<uint64_t, libMesh::Node*> nastran_to_libmesh_node_map;
        /// Mapping from libMesh/MAST node pointers to Nastran grid IDs from BDF input.
        std::map<const libMesh::Node*, uint64_t> libmesh_to_nastran_node_map;
        /// Mapping from Nastran element IDs from BDF input to pointers to libMesh/MAST elements.
        std::map<uint64_t, libMesh::Elem*> nastran_to_libmesh_elem_map;
        /// Mapping from libMesh/MAST element pointers to Nastran element IDs from BDF input.
        std::map<libMesh::Elem*, uint64_t> libmesh_to_nastran_elem_map;

        /**
         * Mapping from Nastran property ID/element-type pair to libMesh/MAST subdomain.
         *
         * Note we must use a pair combining a Nastran property ID along with a specific element
         * type to get a unique map key to libMesh/MAST subdomains. This is because Nastran allows
         * one property ID to be assigned to elements of different type (ie. PSHELL to both CQUAD4
         * and CTRIA3), but libMesh does not allow elements of different type to belong to the
         * same subdomain. In this case, we use multiple subdomains to contain the different
         * element types and store this mapping for reference.
         */
        std::map<std::pair<int, int>, int> nastran_pid_elemtype_to_libmesh_subdomain_map = {};

        void read_nodes(BDFModel* model, libMesh::MeshBase& the_mesh);
        void read_elements(BDFModel* model, libMesh::MeshBase& the_mesh);
        void read_node_boundaries(BDFModel* model, libMesh::MeshBase& the_mesh);

        /**
         * Map from Nastran elements to equivalent libMesh/MAST element types.
         * TODO: Not yet complete, need to add all Nastran elements we need support for.
         */
        std::map<std::string, libMesh::ElemType> nastran_to_libmesh_elem_type_map = {

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

        void initialize_python();
        void finalize_python();
    };

    // Utility functions.
    void printElementMap(std::map<std::string, std::vector<std::vector<int>>> elementMap);
    void printNodeCoords(std::vector<std::vector<double>> nodes);
}

#endif //  MAST_ENABLE_NASTRANIO
#endif // __mast_nastran_io_h__
