#ifndef MAST_SOLID_1D_L_SECTION_ELEMENT_PROPERTY_CARD_H_INCLUDED
#define MAST_SOLID_1D_L_SECTION_ELEMENT_PROPERTY_CARD_H_INCLUDED


#include "property_cards/solid_1d_L_section_element_property_card.h"
#include "material/mast_isotropic_material.h"

extern libMesh::LibMeshInit* p_global_init;

namespace TEST
{
    class AluminumLSection : public MAST::Solid1DLSectionElementPropertyCard
    {
    public:
        
        AluminumLSection(const uint n_target_elems=3500):
        DIM1("DIM1", 6.000), DIM2("DIM2", 3.000), DIM3("DIM3", 0.250),
        DIM4("DIM4", 0.125), offset_y("offy_param", -0.918), 
        offset_z("offz_param", -0.347), DIM1_f("DIM1", DIM1), 
        DIM2_f("DIM2", DIM2), DIM3_f("DIM3", DIM3), DIM4_f("DIM4", DIM4),
        offsety_f("hy_off", offset_y), offsetz_f("hz_off", offset_z),
        centroid_true(xc_true, yc_true), shear_center_true(xs_true, ys_true)
        {
            add(DIM1_f);
            add(DIM2_f);
            add(DIM3_f);
            add(DIM4_f);
            
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
        MAST::Parameter DIM3;
        MAST::Parameter DIM4;
        
        MAST::Parameter offset_y;
        MAST::Parameter offset_z;
        
        MAST::ConstantFieldFunction DIM1_f;
        MAST::ConstantFieldFunction DIM2_f;
        MAST::ConstantFieldFunction DIM3_f;
        MAST::ConstantFieldFunction DIM4_f;
        
        MAST::ConstantFieldFunction offsety_f;
        MAST::ConstantFieldFunction offsetz_f;
        
        // True values
        const Real area_true = 1.8437500000000016e+00;
        const Real torsion_constant_true = 3.2286682945239953e-02;
        const Real first_area_moment_z_true = -1.1769375000000004e+00;
        const Real first_area_moment_y_true = 3.7664687499999974e+00;
        const Real second_area_moment_zz_true = 1.6049689895833241e+00;
        const Real second_area_moment_yy_true = 1.4607873559895829e+01;
        const Real second_area_moment_zy_true = -3.6365401874999983e+00;
        const Real second_area_moment_polar_true = 1.6212842549479152e+01;
        const Real Izzc_true = 8.5368390271891748e-01;
        const Real Iyyc_true = 6.9136162881797452e+00;
        const Real Izyc_true = -1.2322563559322042e+00;
        const Real Ipc_true = 7.7673001908986627e+00;
        const Real warping_constant_true = 9.4267903516538354e-02;
        const Real kappa_z_true = 6.9387360587828650e-01;
        const Real kappa_y_true = 1.5217240110756519e-01;
        const Real xs_true = -3.0394644727391729e-01;
        const Real ys_true = -9.2444644346018534e-01;
        const Real xst_true = -3.0394644727391729e-01;
        const Real yst_true = -9.2444644346018534e-01;
        const Real xc_true = 2.0428305084745730e+00;
        const Real yc_true = -6.3833898305084713e-01;
        const std::vector<libMesh::Point> stress_points_true = {
            libMesh::Point(-1.9803305084745730e+00, 3.5133389830508470e+00, 0.),
            libMesh::Point(3.8946694915254270e+00, 5.1333898305084713e-01, 0.),
            libMesh::Point(-2.1053305084745730e+00, 5.1333898305084713e-01, 0.),
            libMesh::Point(-2.1053305084745730e+00, 3.5133389830508470e+00, 0.)
        };
        
        libMesh::Point centroid_true;
        libMesh::Point shear_center_true;

    };
}
#endif // MAST_SOLID_1D_L_SECTION_ELEMENT_PROPERTY_CARD_H_INCLUDED
