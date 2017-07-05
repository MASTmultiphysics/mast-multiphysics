from mast_interface cimport add_node
from mast_interface cimport add_nothing

cdef public void add_a_lot_of_nodes(void * d):
    add_nothing(d)
    add_node(1, 2, 3, 4, d)

