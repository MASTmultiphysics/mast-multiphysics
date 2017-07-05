cdef extern from "mast_interface.h":
    void add_node(int i, double x, double y, double z, void* d)
    void add_nothing(void* d)
