from pyNastran.bdf.bdf import BDF
from mast_interface cimport Model
#from mast_interface cimport add_node
#from mast_interface cimport add_beam

cdef public void read_and_initialize_model(Model* m):
    myBDF = BDF()
    #nm = '/Users/manav/Documents/codes/numerical_lib/pynastran/pyNastran/models/beam_modes/cbarao_cbeam_static.bdf'
    nm = '/Users/manav/Documents/codes/numerical_lib/pynastran/pyNastran/models/plate_py/plate_py.dat'
    myBDF.read_bdf(nm)

    # iterate over various quantities and add them to the model
    
    ######################################################
    #            NODES
    ######################################################
    # iterate over all the nodes and add them to the model
    for n_id, nd in myBDF.nodes.items():
        xyz = nd.get_position()
        m.add_node(n_id, xyz[0], xyz[1], xyz[2])

    ######################################################
    #            ELEMENTS
    ######################################################
    for e_id, el in myBDF.elements.items():
        type = el.type
        if type == 'CBAR':
            m.add_edge2(e_id, el.Pid(), el.Ga(), el.Gb())
        elif type == 'CQUAD4':
            nds = el.node_ids
            m.add_quad4(e_id, el.Pid(), nds[0], nds[1], nds[2], nds[3])
        else:
            s = "Elem not handled: " + type
            print(s)

    ######################################################
    #            MATERIAL
    ######################################################
    for m_id, mat in myBDF.materials.items():
        type = mat.type
        if type == 'MAT1':
            m.add_isotropic_material(m_id, mat.E(), mat.Nu(), mat.Rho(), mat.a)
        else:
            s = "Material not handled: " + type
            print(s)

    ######################################################
    #            PROPERTY
    ######################################################
    for p_id, p in myBDF.properties.items():
        type = p.type
        if type == 'PSHELL':
            m.add_2d_section_property(p_id, p.Mid1(), p.Thickness())
        else:
            s = "Property not handled: " + type
            print(s)
