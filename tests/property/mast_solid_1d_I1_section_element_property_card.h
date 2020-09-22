#ifndef MAST_SOLID_1D_I1_SECTION_ELEMENT_PROPERTY_CARD_H_INCLUDED
#define MAST_SOLID_1D_I1_SECTION_ELEMENT_PROPERTY_CARD_H_INCLUDED


#include "property_cards/solid_1d_I1_section_element_property_card.h"
#include "material/mast_isotropic_material.h"

extern libMesh::LibMeshInit* p_global_init;

namespace TEST
{
    class AluminumI1Section : public MAST::Solid1DI1SectionElementPropertyCard
    {
    public:
        AluminumI1Section(const uint n_target_elems=3500):
        DIM1("DIM1", 3.770), DIM2("DIM2", 0.170), DIM3("DIM3", 5.470),
        DIM4("DIM4", 5.900), offset_y("offy_param", -0.787), 
        offset_z("offz_param", 0.687), DIM1_f("DIM1", DIM1), 
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
        const Real area_true = 2.6241000000000039e+00;
        const Real torsion_constant_true = 3.4912416307463445e-02;
        const Real first_area_moment_z_true = -2.0651667000000016e+00;
        const Real first_area_moment_y_true = 1.8027567000000015e+00;
        const Real second_area_moment_zz_true = 1.7639240550400139e+01;
        const Real second_area_moment_yy_true = 3.4324069553999919e+00;
        const Real second_area_moment_zy_true = -1.4187695228999926e+00;
        const Real second_area_moment_polar_true = 2.1071647505800129e+01;
        const Real Izzc_true = 1.6013954357500140e+01;
        const Real Iyyc_true = 2.1939131024999914e+00;
        const Real Izyc_true = 7.5495165674510645e-15;
        const Real Ipc_true = 1.8207867460000131e+01;
        const Real warping_constant_true = 1.7690289289249652e+01;
        const Real kappa_z_true = 5.4403106725573769e-01;
        const Real kappa_y_true = 3.5065087923352467e-01;
        const Real xs_true = 6.8699874537432160e-01;
        const Real ys_true = -7.8699671905199586e-01;
        const Real xst_true = 6.8699874537432160e-01;
        const Real yst_true = -7.8699671905199586e-01;
        const Real xc_true = 6.8699999999999961e-01;
        const Real yc_true = -7.8699999999999948e-01;
        const std::vector<libMesh::Point> stress_points_true = {
            libMesh::Point(1.9700000000000000e+00, 2.9500000000000002e+00, 0.),
            libMesh::Point(1.9700000000000000e+00, -2.9500000000000002e+00, 0.),
            libMesh::Point(-1.9700000000000000e+00, -2.9500000000000002e+00, 0.),
            libMesh::Point(-1.9700000000000000e+00, 2.9500000000000002e+00, 0.)
        };
        
        libMesh::Point centroid_true;
        libMesh::Point shear_center_true;

    };
}
#endif // MAST_SOLID_1D_I1_SECTION_ELEMENT_PROPERTY_CARD_H_INCLUDED
