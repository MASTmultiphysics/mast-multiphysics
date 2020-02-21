#ifndef __mast__cross_section__
#define __mast__cross_section__

// libMesh Includes
#include "libmesh/libmesh.h"
#include "libmesh/replicated_mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/point.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/equation_systems.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/mesh_triangle_interface.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/linear_solver.h"

// // MAST includes
#include "base/constant_field_function.h"
#include "property_cards/solid_1d_section_element_property_card.h"
#include "base/warping_assembly.h"
#include "elasticity/warping_system_initialization.h"
#include "elasticity/structural_nonlinear_assembly.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "property_cards/solid_2d_section_element_property_card.h"
#include "base/physics_discipline_base.h"



namespace MAST 
{
    
class CrossSection
{
    public:
        CrossSection(const libMesh::LibMeshInit& init, 
                     const uint target_n_elems,
                     MAST::Solid1DSectionElementPropertyCard& property,
                     const libMesh::ElemType element_type=libMesh::TRI6);
        
        /*!
        *   virtual destructor
        */
        virtual ~CrossSection() { }
        
        void calculate_shoelace_geometric_properties();
        void calculate_geometric_properties();
        void calculate_warping_properties();
        
        Real get_area(const libMesh::Point& p, const Real t);
        Real get_area_derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t);
        
        libMesh::Point get_centroid(const libMesh::Point& p, const Real t);
        libMesh::Point get_centroid_derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t);
        
        Real        get_first_area_moment_z(const libMesh::Point& p, const Real t);
        Real        get_first_area_moment_z_derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t);
        
        Real        get_first_area_moment_y(const libMesh::Point& p, const Real t);
        Real        get_first_area_moment_y_derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t);
        
        RealVectorX get_second_area_moments(const libMesh::Point& p, const Real t);
        RealVectorX get_second_area_moments_derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t);
        
        Real        get_second_polar_area_moment(const libMesh::Point& p, const Real t);
        Real        get_second_polar_area_moment_derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t);
                
        Real get_torsion_constant(const libMesh::Point& p, const Real t);
        Real get_torsion_constant_derivative(MAST::FunctionBase& f, const libMesh::Point& p, const Real t);
        
        libMesh::Point get_shear_center(const libMesh::Point& p, const Real t);
        libMesh::Point get_shear_center_derivative(MAST::FunctionBase& f, const libMesh::Point& p, const Real t);
        
        RealVectorX get_shear_coefficients(const libMesh::Point& p, const Real t);
        RealVectorX get_shear_coefficients_derivative(MAST::FunctionBase& f, const libMesh::Point& p, const Real t);
        
        Real get_warping_constant(const libMesh::Point& p, const Real t);
        Real get_warping_constant_derivative(MAST::FunctionBase& f, const libMesh::Point& p, const Real t);
                
        RealVectorX get_second_area_moments_about_centroid(const libMesh::Point& p, const Real t);
        RealVectorX get_radii_of_gyration(const libMesh::Point& p, const Real t);
        RealVectorX get_elastic_section_modulii(const libMesh::Point& p, const Real t);
        

    protected:
        
        /*!
         *  METHODS
         */
        
        void _calculate_shoelace_geometric_properties();
        
        void _calculate_shoelace_geometric_properties_sensitivities(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t);
        
        bool _shoelace_calculated = false;
        
        void _configure_linear_solver();
        
        void _mesh_cross_section(const libMesh::Point& p, const Real t);
        
        void _build_model(const libMesh::Point& p, const Real t);
        
        void _solve_warping();
        
        void _update_geometry(const libMesh::Point& p, const Real t);
        
        void _update_geometry(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t);
        
        void _update_warping_geometry(const libMesh::Point& p, const Real t);
        
        
        /*!
         * VARIABLES
         */
        const libMesh::LibMeshInit& _init;
        const uint _target_n_elems;
        const libMesh::ElemType _element_type;
        
        std::unique_ptr<libMesh::ReplicatedMesh> _mesh;
        std::unique_ptr<libMesh::EquationSystems> _equation_systems;
        std::unique_ptr<libMesh::FEType> _fetype;
        MAST::NonlinearSystem* _system;
        std::unique_ptr<MAST::PhysicsDisciplineBase> _discipline;
        std::unique_ptr<MAST::WarpingSystemInitialization> _warping_system;
        std::unique_ptr<MAST::WarpingAssembly> _assembly;
        std::unique_ptr<MAST::StructuralNonlinearAssemblyElemOperations> _elem_ops;
        MAST::NonlinearSystem* _nonlinear_system;
        
        // Define geometric properties as MAST Parmeters
        std::map<std::string, std::unique_ptr<MAST::Parameter>> _parameters;
        std::map<std::string, std::unique_ptr<MAST::ConstantFieldFunction>> _field_functions;
        
        // Create the section property card and attach all geometric property values.
        std::unique_ptr<MAST::Solid2DSectionElementPropertyCard> _section;
        
        MAST::Solid1DSectionElementPropertyCard& _property;
        
        Real _A, _Cx, _Cy, _Qx, _Qy, _Ixx, _Iyy, _Ixy, _rx, _ry, _Ixxc, _Iyyc, 
             _Ixyc, _Z_xx_p, _Z_xx_n, _Z_yy_p, _Z_yy_n;
             
        Real _dA, _dCx, _dCy, _dQx, _dQy, _dIxx, _dIyy, _dIxy, _drx, _dry, 
             _dIxxc, _dIyyc, _dIxyc, _dZ_xx_p, _dZ_xx_n, _dZ_yy_p, _dZ_yy_n;
        
        std::unique_ptr<libMesh::NumericVector<double>> _F_warp;
        std::unique_ptr<libMesh::NumericVector<double>> _F_shearx;
        std::unique_ptr<libMesh::NumericVector<double>> _F_sheary;
        
        std::unique_ptr<libMesh::NumericVector<Real>> _Omega;
        std::unique_ptr<libMesh::NumericVector<Real>> _Psi;
        std::unique_ptr<libMesh::NumericVector<Real>> _Phi;
        
        std::vector<libMesh::Point> _geometry_points;
        std::vector<std::vector<libMesh::Point>> _hole_list;
        std::string _sensitivity_parameter;
    //     std::vector<libMesh::Point> _stress_points;
        
        bool _warping_calculated = false;
        
        geometric_properties _gp;
        warping_properties _wp;
};
    
}

#endif // __mast__cross_section__
