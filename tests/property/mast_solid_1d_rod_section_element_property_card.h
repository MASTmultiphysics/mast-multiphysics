#ifndef MAST_SOLID_1D_ROD_SECTION_ELEMENT_PROPERTY_CARD_H_INCLUDED
#define MAST_SOLID_1D_ROD_SECTION_ELEMENT_PROPERTY_CARD_H_INCLUDED

#include "property_cards/solid_1d_rod_section_element_property_card.h"
#include "material/mast_isotropic_material.h"

extern libMesh::LibMeshInit* p_global_init;

#define PI 3.1415926535897932

namespace TEST
{
    class AluminumRodSection : public MAST::Solid1DRodSectionElementPropertyCard
    {
    public:
        
        AluminumRodSection(const uint n_target_elems=3500):
        DIM1("DIM1", 3.234),
        offset_y("offy_param", 0.287), offset_z("offz_param", -1.654), 
        DIM1_f("DIM1", DIM1),
        offsety_f("hy_off", offset_y), offsetz_f("hz_off", offset_z),
        centroid_true(xc_true, yc_true), shear_center_true(xs_true, ys_true)
        {
            add(DIM1_f);
            
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
        
        MAST::Parameter offset_y;
        MAST::Parameter offset_z;
        
        MAST::ConstantFieldFunction DIM1_f;
        
        MAST::ConstantFieldFunction offsety_f;
        MAST::ConstantFieldFunction offsetz_f;
        
        // True values calculated using the sectionproperties module in python BAR.py
        const Real r = DIM1();
        const Real area_true = PI*r*r;
        const Real torsion_constant_true = PI*pow(r, 4.0)/2.0;
        const Real first_area_moment_z_true = area_true * offset_y();
        const Real first_area_moment_y_true = area_true * offset_z();
        const Real Izzc_true = PI*pow(r, 4.0)/4.0;
        const Real Iyyc_true = PI*pow(r, 4.0)/4.0;
        const Real Izyc_true = 0.0;
        const Real Ipc_true = PI*pow(r, 4.0)/2.0;
        const Real second_area_moment_zz_true = Izzc_true + area_true * pow(offset_y(), 2.);
        const Real second_area_moment_yy_true = Iyyc_true + area_true * pow(offset_z(), 2.);
        const Real second_area_moment_zy_true = Izyc_true + area_true * offset_y() * offset_z();
        const Real second_area_moment_polar_true = second_area_moment_zz_true + second_area_moment_yy_true;
        const Real warping_constant_true = 0.0;
        const Real kappa_z_true = 8.4967018479575718e-01;
        const Real kappa_y_true = 8.4967018474199052e-01;
        const Real xs_true = offset_z();
        const Real ys_true = offset_y();
        const Real xc_true = offset_z();
        const Real yc_true = offset_y();
        const std::vector<libMesh::Point> stress_points_true = {
            libMesh::Point(0., r, 0.),
            libMesh::Point(r, 0., 0.),
            libMesh::Point(0., -r, 0.),
            libMesh::Point(-r, 0., 0.)
        };

        
        libMesh::Point centroid_true;
        libMesh::Point shear_center_true;

    };
}
#endif // MAST_SOLID_1D_ROD_SECTION_ELEMENT_PROPERTY_CARD_H_INCLUDED
