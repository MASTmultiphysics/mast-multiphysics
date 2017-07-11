from pyNastran.bdf.bdf import BDF
from mast_interface cimport Model
from libcpp.string cimport string




cdef public void read_and_initialize_model(Model* m, string f_nm):
    myBDF = BDF()
    nm = f_nm.decode("UTF-8")
    myBDF.read_bdf(nm)

    #tell the solver about the kind of solution specified in the file
    m.set_sol(myBDF.sol)

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
        t = el.type
        if t == 'CBAR':
            m.add_edge2(e_id, el.Pid(), el.Ga(), el.Gb())
        elif t == 'CQUAD4':
            nds = el.node_ids
            m.add_quad4(e_id, el.Pid(), nds[0], nds[1], nds[2], nds[3])
        else:
            s = "Elem not handled: " + type
            print(s)
            assert False


    ######################################################
    #            MATERIAL
    ######################################################
    for m_id, mat in myBDF.materials.items():
        t = mat.type
        if t == 'MAT1':
            m.add_isotropic_material(m_id, mat.E(), mat.Nu(), mat.Rho(), mat.a)
        else:
            s = "Material not handled: " + type
            print(s)
            assert False

    ######################################################
    #            PROPERTY
    ######################################################
    for p_id, p in myBDF.properties.items():
        t = p.type
        if t == 'PSHELL':
            m.add_2d_section_property(p_id, p.Mid1(), p.Thickness())
        else:
            s = "Property not handled: " + type
            print(s)
            assert False

    ######################################################
    #            SubCases
    ######################################################
    for s_id, s in myBDF.subcases.items():
        if s_id > 0:  # PyNastran seems to be reading the header as subcase 0
            load  = s.params['LOAD'][0]   #load set ids
            spc   = s.params['SPC'][0]    #spc set ids
            title = s.params['TITLE'][0]
            m.add_subcase(s_id, load, spc);






cdef void add_spc_to_model(spc, Model* m):
    t = spc.type
    if t == 'SPC':
        i = 0
        for n in spc.node_ids:
            dofs = spc.components[i]
            val  = spc.enforced[i]
            for c in dofs:
                m.add_spc(n, int(c), val)
            i = i+1
    elif t == 'SPC1':
        nodes = spc.node_ids
        dofs  = spc.constraints
        for c in dofs:
            # tell the model to constraint this dof for the nodes
            for n in nodes:
                m.add_spc(n, int(c), 0.)
    else:
        str = "SPC Type not handled: " + type
        print(str)
        assert False





cdef void process_spc_list(spc, Model* m):
    # make sure that the spcs are a list
    assert isinstance(spc, list)
    for s in spc:
        t = s.type
        if t == 'SPC' or t == 'SPC1':
            add_spc_to_model(s, m)
        elif t == 'SPCADD':
            # this will iterate on all the BCs in the i_th set. pyNastran
            # uses a list to store all BCs in a BC set
            for i in s.sets:
                process_spc_list(i, m)
        else:
            str = "SPC Type not handled: " + type
            print(str)
            assert False





cdef void add_load_to_model(scale, load, Model* m):
    t = load.type
    if t == 'FORCE':
        ff = load.transform_load()
        #assert ff[0]  #True if the magnitude is greater than 0.
        n   = ff[1]  # node id
        f   = ff[2]  # force vector in global coordinate system
        m.add_force(n, scale*f[0], scale*f[1], scale*f[2])
    else:
        str = "Load Type not handled: " + type
        print(str)
        assert False






cdef void process_load_list(scale, loads, Model* m):
    # make sure that the loads are a list
    assert isinstance(loads, list)
    for l in loads:
        t = l.type
        if t == 'FORCE':
            add_load_to_model(scale, l, m)
        elif t == 'LOAD':
            l_scale = l.scale  # the scaling factor of the whole load card
            i = 0
            # this will iterate on all the loads in the i_th set. pyNastran
            # uses a list to store all loads in a load set
            for lset in l.load_ids:
                # combine the load factor for this load case with that of the
                # whole load card
                process_load_list(scale*l_scale*l.scale_factors[i], lset, m)
                i = i+1
        else:
            str = "Load Type not handled: " + type
            print(str)
            assert False






cdef public void read_and_initialize_loads_for_subcase(int subid, Model* m, string f_nm):
    myBDF    = BDF()
    nm       = f_nm.decode("UTF-8")
    myBDF.read_bdf(nm)
    subcase  = myBDF.subcases[subid]
    load_id  = subcase.params['LOAD'][0]   #load set ids
    spc_id   = subcase.params['SPC'][0]    #spc set ids
    
    
    ######################################################
    #            SPC
    ######################################################
    spc = myBDF.spcs[spc_id]
    process_spc_list(spc, m)

    ######################################################
    #            LOADS
    ######################################################
    loads = myBDF.loads[load_id]
    process_load_list(1., loads, m) # use a unity scaling



