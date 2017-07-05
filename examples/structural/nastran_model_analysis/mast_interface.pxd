cdef extern from "examples/structural/nastran_model_analysis/mast_interface.h" namespace "MAST":
    cdef cppclass Model:
        void add_node(int n_id, double x, double y, double z)
        void add_edge2(int e_id, int p_id,   int n1,   int n2)
        void add_quad4(int e_id, int p_id,   int n1,   int n2,   int n3,   int n4)
        void add_isotropic_material(int m_id, double E, double nu, double rho, double alpha)
        void add_2d_section_property(int p_id, int m_id, double t)
