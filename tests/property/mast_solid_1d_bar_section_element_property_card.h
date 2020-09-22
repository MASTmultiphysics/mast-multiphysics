#ifndef MAST_SOLID_1D_BAR_SECTION_ELEMENT_PROPERTY_CARD_H_INCLUDED
#define MAST_SOLID_1D_BAR_SECTION_ELEMENT_PROPERTY_CARD_H_INCLUDED

#include "property_cards/solid_1d_bar_section_element_property_card.h"
#include "material/mast_isotropic_material.h"

extern libMesh::LibMeshInit* p_global_init;

namespace TEST
{
    class AluminumBarSection : public MAST::Solid1DBarSectionElementPropertyCard
    {
    public:
        
        AluminumBarSection(const uint n_target_elems=3500):
        DIM1("DIM1", 3.0), DIM2("DIM2", 0.75), 
        offset_y("offy_param", 0.287), offset_z("offz_param", 1.654), 
        DIM1_f("DIM1", DIM1), DIM2_f("DIM2", DIM2),
        offsety_f("hy_off", offset_y), offsetz_f("hz_off", offset_z),
        centroid_true(xc_true, yc_true), shear_center_true(xs_true, ys_true)
        {
            add(DIM1_f);
            add(DIM2_f);
            
            add(offsety_f);
            add(offsetz_f);
            
            // Add the material card to the section card
            set_material(material);
            
            // Specify a section orientation point and add it to the 
            RealVectorX orientation = RealVectorX::Zero(3);
            orientation(1) = 1.0;
            y_vector() = orientation;
            
            // Now initialize the section
            init(*p_global_init, n_target_elems);
        }
        
        TEST::Aluminum7075T6 material;
        
        MAST::Parameter DIM1;
        MAST::Parameter DIM2;
        
        MAST::Parameter offset_y;
        MAST::Parameter offset_z;
        
        MAST::ConstantFieldFunction DIM1_f;
        MAST::ConstantFieldFunction DIM2_f;
        
        MAST::ConstantFieldFunction offsety_f;
        MAST::ConstantFieldFunction offsetz_f;
        
        // True values calculated using the sectionproperties module in python BAR.py
        const Real area_true = 2.2500000000000004e+00;
        const Real torsion_constant_true = 3.5540414988249380e-01;
        const Real first_area_moment_z_true = 6.4574999999999994e-01;
        const Real first_area_moment_y_true = 3.7214999999999945e+00;
        const Real second_area_moment_zz_true = 2.9079899999999986e-01;
        const Real second_area_moment_yy_true = 7.8428609999999814e+00;
        const Real second_area_moment_zy_true = 1.0680705000000006e+00;
        const Real second_area_moment_polar_true = 8.1336599999999812e+00;
        const Real Izzc_true = 1.0546874999999994e-01;
        const Real Iyyc_true = 1.6875000000000009e+00;
        const Real Izyc_true = 2.4424906541753444e-15;
        const Real warping_constant_true = 6.1030611538894504e-02;
        const Real kappa_z_true = 8.3330248143102326e-01;
        const Real kappa_y_true = 5.5763491072129889e-01;
        const Real xs_true = 1.6540000458463804e+00;
        const Real ys_true = 2.8700000023051281e-01;
        const Real xst_true = 1.6540000458463804e+00;
        const Real yst_true = 2.8700000023051281e-01;
        const Real xc_true = 1.654;
        const Real yc_true = 0.287;
        const std::vector<libMesh::Point> stress_points_true = {
            libMesh::Point(1.5,     0.375, 0.),
            libMesh::Point(1.5,     -0.375, 0.),
            libMesh::Point(-1.5,    -0.375, 0.),
            libMesh::Point(-1.5,    0.375, 0.)
        };
        
        libMesh::Point centroid_true;
        libMesh::Point shear_center_true;

    };
}
#endif // MAST_SOLID_1D_BAR_SECTION_ELEMENT_PROPERTY_CARD_H_INCLUDED
