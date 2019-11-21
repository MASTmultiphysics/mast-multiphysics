# distutils: language = c++
# ^^^ This line above is read by setup.py and lets it know to compile the
# program as C++ instead of the default, C.

"""
pyNASTRAN Notes
---------------
Dependencies: cpylog, scipy, numpy
"""

# pyNastran Imports
from pyNastran.bdf.bdf import BDF

# Standard C++ Imports
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.map cimport map
from libcpp cimport bool
from libcpp.list cimport list as clist
from libcpp.set cimport set as cset

# Cython Imports
from cython.operator import dereference, postincrement
from warnings import warn

# Setup Numpy Usage
import numpy as np
cimport numpy as np
DTYPE = np.float
ctypedef np.float_t DTYPE_t


cdef public class BDFModel[object BDFModel, type BDFModelType]:
    cdef:
        string bdfPath

    # The readonly specifier allows C/C++ code to still access these variables, but prevents them from changing it.
    cdef readonly :
        int nDims
        int nNodes
        int nElems
        int nMaterials
        int nProperties
        map[string, int] nElemTypes

    # Need to declare this here, even without a type specifier. Otherwise,
    # a "Segmentation fault (core dumped)" error occurs.
    cdef myBDF

    def __cinit__(self, string bdf_path):
        self.bdfPath = bdf_path
        self.myBDF = BDF()
        nm = bdf_path.decode("UTF-8")
        self.myBDF.read_bdf(nm)
        self.nNodes = self.myBDF.nnodes
        self.nElems = self.myBDF.nelements
        self.nMaterials = self.myBDF.nmaterials
        self.nProperties = self.myBDF.nproperties
        # self.temperatures = None

        # Perturb the nodes of zero-length vectors, this is a hack to make zero-length
        # elements work in MAST.
        perturbZeroLengthBushings(self)


cdef public BDFModel buildBDFModel(string bdf_path):
    """
    Acts as the constructor for the BDFModel class.
    """
    return BDFModel(bdf_path)


cdef public void printBDFStats(BDFModel model):
    print(model.myBDF.get_bdf_stats())
    print("card_count   = %s\n"%model.myBDF.card_count)
    print("reject_count = %s\n"%model.myBDF.reject_count)


cdef public vector[vector[double]] getNodes(BDFModel model):
    """
    Creates a vector of vectors where each row is the node id and coordinates of
    a single node (id, x, y, z). Also determines the dimension of the model.
    
    Notes
    ----
    Have to manually add "#include <vector>" to the pynastran.h file for this 
    to work. Is there a way to tell setup.py to automatically do this?
    
    References
    ----------
    "Cython, confusion, regarding the vector (size) constructor", StackOverflow,
    Ami Tavory, Jan 16, 2016.
    """
    # Get the number of rows and columns
    cdef int n = model.nNodes # Number of Nodes
    cdef int m = 4            # node number + 3 dimensions = 4

    # Create a 2D vector (vector of vectors) to store node numbers and coords
    cdef vector[double] node_vec
    cdef vector[vector[double]] nodes
    node_vec = vector[double](m)
    nodes = vector[vector[double]](n, node_vec)

    # Iterate over all nodes and all them to the 2D vector
    cdef bool anyZ = False
    cdef bool anyY = False
    cdef bool anyX = False
    cdef int i = 0
    for key in model.myBDF.nodes:
        node = model.myBDF.nodes[key]
        nodes[i][0] = node.nid
        nodes[i][1] = node.xyz[0]
        nodes[i][2] = node.xyz[1]
        nodes[i][3] = node.xyz[2]
        if node.xyz[0]!=0.0:
            anyX = True
        if node.xyz[1]!=0.0:
            anyY = True
        if node.xyz[2]!=0.0:
            anyZ = True
        i += 1

    if anyZ:
        model.nDims = 3
    elif anyY:
        model.nDims = 2
    elif anyX:
        model.nDims = 1
    else:
        model.nDims = 0

    return nodes


cdef public map[string, vector[vector[int]]] getElements(BDFModel model):
    """
    Creates a map of the elements contained in the model. The keys of the map are the element type and the values
    of the map is a 2D vector (vector of vector<int>) which stores element ID, property ID, and node numbers.
    """
    cdef int n = model.nElems     # Number of elements

    cdef map[string, vector[vector[int]]] elementMap;

    # Pre-allocate space to store the element information
    cdef map[string, int] elemTypes = getNumberOfElementTypes(model)
    cdef map[string, int] i_counter
    cdef map[string, int].iterator it = elemTypes.begin()
    cdef vector[int] elem_vec
    cdef vector[vector[int]] elements
    cdef int m, nElemNodes
    while(it != elemTypes.end()):
        # print(dereference(it).first) # print the key
        # print(dereference(it).second) # print the associated value
        elemType = dereference(it).first
        nElemTypes = dereference(it).second
        nElemNodes = int(str(elemType).split("_")[-1].strip("'"))
        m = nElemNodes+2

        elem_vec = vector[int](m)
        elementMap[elemType] = vector[vector[int]](nElemTypes, elem_vec)

        i_counter[elemType] = 0

        postincrement(it) # Increment the iterator to the net element
    mapKeys = []

    # Loop through the elements in the BDF and all the to element map
    cdef int i
    cdef int j
    for key,element in model.myBDF.elements.items():
        if "offset" in element.object_attributes():
            if element.offset!=0:
                raise NotImplementedError("Support of non-zero offsets not yet implemented!")
        if "zoffset" in element.object_attributes():
            if element.zoffset!=0:
                raise NotImplementedError("Support of non-zero offets not yet implemented!")
        nElemNodes = len(element.nodes)
        elemType = (element.type + "_%i"%(nElemNodes)).encode('utf-8')
        i = i_counter[elemType]
        elementMap[elemType][i][0] = element.eid
        elementMap[elemType][i][1] = element.pid
        for j in range(nElemNodes):
            elementMap[elemType][i][j+2] = element.nodes[j]
        i_counter[elemType] = i_counter[elemType]+1

    return elementMap


cdef public map[string, int] getNumberOfElementTypes(BDFModel model):
    """
    Gets a map of the different element types and the number of each element type in the model. Useful for 
    pre-allocation when importing elements.
    """
    cdef map[string, int] elementTypes
    elemTypeList = []
    for key, element in model.myBDF.elements.items():
        elemType = (element.type + "_%i"%(len(element.nodes))).encode("utf-8")
        if elemType not in elemTypeList:
            elemTypeList.append(elemType)
            elementTypes[elemType] = 1
        else:
            elementTypes[elemType] += 1

    model.nElemTypes = elementTypes
    return elementTypes


cdef perturbZeroLengthBushings(BDFModel model):
    """
    As of November 18, 2019, MAST does not support the specification of the elemental x-axis which is required when 
    the element has a zero-length (i.e. CBUSH in Nastran) and the x-axis cannot be determined from the geometry. This
    method provides a work around by getting the x-axis from the coordinate system the bushing references and then
    perturbing one of the nodes by a very small distance in the direction of the x-axis. This then allows MAST to 
    determine the x-axis from geometry and introduces a negligible effect on the stiffness matrix of the bushing
    element.
    """
    cdef int k
    cdef double eps, nzero, L
    eps=1.4901161193847656e-08
    nzero=1e-08

    # Get element IDs of bushing elements
    data = model.myBDF.get_elements_properties_nodes_by_element_type()
    eids = []
    if "CBUSH" in data.keys():
        eids += (data["CBUSH"][0]).tolist()
    if "CBUSH1D" in data.keys():
        eids += (data["CBUSH1D"][0]).tolist()
    if "CBUSH2D" in data.keys():
        eids += (data["CBUSH2D"][0]).tolist()

    # Loop through all bushing elements
    for k in range(len(eids)):
        eid = eids[k]
        bushing = model.myBDF.elements[eid]
        L = np.linalg.norm(bushing.nodes_ref[1].xyz - bushing.nodes_ref[0].xyz)
        if abs(L)<nzero:
            # Slightly perturb the second node so that L>0 and an x-axis can be defined by the geometry
            i = bushing.cid_ref.i
            di = eps*i
            absmindi = abs(di).min()
            if absmindi<eps: # Imposing a minimum perurbation size to coordinates
                di *= (eps/absmindi)
            model.myBDF.nodes[bushing.nodes[1]].xyz += di
            print("Perturbing GRID %i by %s to define x-axis of zero-length element."%(bushing.nodes[1], di))

            # Define the orientation of the bushing
            bushing.x = bushing.cid_ref.j.tolist()