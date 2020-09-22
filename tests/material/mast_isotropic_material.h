#ifndef MAST_ISOTROPIC_MATERIAL_H_INCLUDED
#define MAST_ISOTROPIC_MATERIAL_H_INCLUDED

#include "property_cards/isotropic_material_property_card.h"


namespace TEST
{    
    class Aluminum7075T6 : public MAST::IsotropicMaterialPropertyCard
    {
    public:
        Aluminum7075T6():
        E("E_param", 71.7e9), nu("nu_param", 0.33), rho("rho_param", 2810.0),
        alpha("alpha_param", 2.52e-05), cp("cp_param", 960.0), 
        k("k_param", 130.0), E_f("E", E), nu_f("nu", nu), rho_f("rho", rho),
        alpha_f("alpha_expansion", alpha), cp_f("cp", cp), k_f("k_th", k)
        {
            add(E_f);
            add(nu_f);
            add(rho_f);
            add(alpha_f);
            add(cp_f);
            add(k_f);
        }
        
        Real G()
        {
            return E() / (2.0 * (1.0 + nu()));
        };
        
        MAST::Parameter E;
        MAST::Parameter nu;
        MAST::Parameter rho;
        MAST::Parameter alpha;
        MAST::Parameter cp;
        MAST::Parameter k;
        
        MAST::ConstantFieldFunction E_f;
        MAST::ConstantFieldFunction nu_f;
        MAST::ConstantFieldFunction rho_f;
        MAST::ConstantFieldFunction alpha_f;
        MAST::ConstantFieldFunction cp_f;
        MAST::ConstantFieldFunction k_f;
    };
    
    
    class Titanium6Al4V : public MAST::IsotropicMaterialPropertyCard
    {
    public:
        Titanium6Al4V():
        E("E_param", 113.8e9), nu("nu_param", 0.342), rho("rho_param", 4430),
        alpha("alpha_param", 9.7e-06), cp("cp_param", 526.3), 
        k("k_param", 6.7), E_f("E", E), nu_f("nu", nu), rho_f("rho", rho),
        alpha_f("alpha_expansion", alpha), cp_f("cp", cp), k_f("k_th", k)
        {
            add(E_f);
            add(nu_f);
            add(rho_f);
            add(alpha_f);
            add(cp_f);
            add(k_f);
        }
        
        Real G()
        {
            return E() / (2.0 * (1.0 + nu()));
        };
        
        MAST::Parameter E;
        MAST::Parameter nu;
        MAST::Parameter rho;
        MAST::Parameter alpha;
        MAST::Parameter cp;
        MAST::Parameter k;
        
        MAST::ConstantFieldFunction E_f;
        MAST::ConstantFieldFunction nu_f;
        MAST::ConstantFieldFunction rho_f;
        MAST::ConstantFieldFunction alpha_f;
        MAST::ConstantFieldFunction cp_f;
        MAST::ConstantFieldFunction k_f;
    };
}
#endif // MAST_ISOTROPIC_MATERIAL_H_INCLUDED
