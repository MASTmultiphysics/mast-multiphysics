
// C++ includes
#include <iostream>


// MAST includes
#include "heat_conduction/heat_conduction_discipline.h"
#include "heat_conduction/heat_conduction_system_initialization.h"
#include "heat_conduction/heat_conduction_nonlinear_assembly.h"
#include "base/physics_discipline_base.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "boundary_condition/dirichlet_boundary_condition.h"
#include "property_cards/isotropic_material_property_card.h"
#include "property_cards/solid_2d_section_element_property_card.h"
#include "driver/driver_base.h"


// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/numeric_vector.h"



/*!
 *   This defines the property as a linear function of temperature
 *   \f[  p   =   k0 + k1 * temp  \f]
 */
class
ConductanceTempDep:
public MAST::FieldFunction<Real> {
public:
    ConductanceTempDep(const std::string& nm,
                       const MAST::Parameter& k0,
                       const MAST::Parameter& k1,
                       const MAST::FieldFunction<RealVectorX>& temp):
    MAST::FieldFunction<Real>(nm),
    _k0(k0),
    _k1(k1),
    _temp(temp) {
    
    }
    
    ConductanceTempDep(const ConductanceTempDep& f):
    MAST::FieldFunction<Real>(f),
    _k0(f._k0),
    _k1(f._k1),
    _temp(f._temp) {
    
    }
    
    virtual ~ConductanceTempDep() {
        
    }
    

    /*!
     *   @returns a clone of the function
     */
    virtual std::auto_ptr<MAST::FieldFunction<Real> > clone() const {
        MAST::FieldFunction<Real>* rval =
        new ConductanceTempDep(*this);
        
        return std::auto_ptr<MAST::FieldFunction<Real> >(rval);
    }
    
    
    /*!
     *    calculates the value of the function at the specified point,
     *    \par p, and time, \par t, and returns it in \p v.
     */
    virtual void operator() (const libMesh::Point& p,
                             const Real t,
                             Real& v) const {
        Real
        k0   = _k0(),
        k1   = _k1();

        RealVectorX temp(1);
        
        // get temperature value
        _temp(p, t, temp);
        
        // calculate the property value using this temperature
        v = k0 + k1*temp(0);
    }
    
    
    /*!
     *    calculates the value of the function at the specified point,
     *    \par p, and time, \par t, and returns it in \p v.
     */
    virtual void derivative (const MAST::DerivativeType d,
                             const MAST::FunctionBase& f,
                             const libMesh::Point& p,
                             const Real t,
                             Real& v) const {
        Real
        k1   = _k1();
        
        RealVectorX dtemp(1);
        
        // get temperature value
        _temp.derivative(d, f, p, t, dtemp);
        
        // calculate the property sensitivity using this temperature
        v = k1*dtemp(0);
    }
    
    
protected:
    
    const MAST::Parameter& _k0;
    
    const MAST::Parameter& _k1;
    
    const MAST::FieldFunction<RealVectorX>& _temp;
};





int main(int argc, const char * argv[]) {
    
    // initialize the libMesh object
    libMesh::LibMeshInit              init(argc, argv);
    libMesh::ParallelMesh             mesh(init.comm());
    libMesh::EquationSystems          eq_sys(mesh);
    
    // add the system to be used for analysis
    libMesh::System& sys = eq_sys.add_system<libMesh::NonlinearImplicitSystem>("thermal");
    
    // initialize the mesh
    unsigned int
    nx        = 10,
    ny        = 10;
    const Real
    width     = 10.,
    height    = 10.;
    libMesh::MeshTools::Generation::build_square(mesh,
                                                 nx, ny,
                                                 0., width,
                                                 0., height);
    mesh.prepare_for_use();
    
    // variable type
    libMesh::FEType fe_type(libMesh::FIRST,
                            libMesh::LAGRANGE);
    
    MAST::HeatConductionDiscipline            heat_cond(eq_sys);
    MAST::HeatConductionSystemInitialization  heat_cond_sys(sys,
                                                            sys.name(),
                                                            fe_type);
    
    
    // add the Dirichlet boundary condition on left boundary
    MAST::DirichletBoundaryCondition   dirichlet;
    dirichlet.init(3, heat_cond_sys.vars());
    
    // add the boundary condition to the physics
    heat_cond.add_dirichlet_bc(3, dirichlet);
    
    // tell the system about the constraints. This needs to happen
    // before the equation system initialization
    heat_cond.init_system_dirichlet_bc(sys);
    
    // initialize the equation system for analysis
    eq_sys.init();
    
    // print the information
    mesh.print_info();
    eq_sys.print_info();
    
    // add a flux load on the right boundary
    // this parameter defines the constant value of the flux
    MAST::Parameter q ("q", 5.0);
    
    // define the field function for boudnary condition
    MAST::ConstantFieldFunction q_flux("heat_flux", q);
    
    // create a flux boundary condition based on this
    MAST::BoundaryConditionBase flux_load(MAST::HEAT_FLUX);
    flux_load.add(q_flux);
    
    // tell the physics about this load
    heat_cond.add_side_load(0, flux_load);
    
    // initialize the material properties. First the parameters, which
    // are used to define the properties to be constants
    MAST::Parameter k0  (  "k0", 2.0);    // thermal conductance constant coeff
    MAST::Parameter k1  (  "k1", 10.0);   // thermal conductance linear coeff
    
    // Temperature function used to calculate the material property
    MAST::FieldFunction<RealVectorX>& temp = heat_cond_sys.get_solution_function();
    
    // define material property based on the parameter
    ConductanceTempDep  k_th("k_th", k0, k1, temp);
    
    
    // add them to a material property card
    MAST::IsotropicMaterialPropertyCard       materials(1);
    materials.add(k_th);
    
    
    // initialize the element section property
    MAST::Parameter h("h", 0.02);
    
    // define a constant field function for thickness
    MAST::ConstantFieldFunction  h_val("h", h);
    
    // add to a section property card
    MAST::Solid2DSectionElementPropertyCard section_property(1);
    section_property.add(h_val);
    
    // tell the section property to use the material property card
    section_property.set_material(materials);
    
    // tell the conduction physics about the section properties
    heat_cond.set_property_for_subdomain(0, section_property);
    
    // create the nonlinear assembly object
    MAST::HeatConductionNonlinearAssembly   assembly;
    
    // now solve the system
    MAST::Driver::nonlinear_solution(heat_cond, heat_cond_sys, assembly);
    
    // write the solution for visualization
    libMesh::ExodusII_IO(mesh).write_equation_systems("mesh.exo", eq_sys);
    
    return 0;
}



//
//// C++ includes
//#include <iostream>
//
//
//// MAST includes
//
//
//// libMesh includes
//#include "libmesh/libmesh.h"
//#include "libmesh/parallel_mesh.h"
//#include "libmesh/abaqus_io.h"
//#include "libmesh/exodusII_io.h"
//
//
//int main(int argc, const char * argv[]) {
//    
//    // initialize the libMesh object
//    libMesh::LibMeshInit              init(argc, argv);
//    libMesh::ParallelMesh             mesh(init.comm());
//    
//    // initialize the mesh
//    libMesh::AbaqusIO abaqus(mesh);
//    abaqus.build_sidesets_from_nodesets = true;
//    abaqus.read("/Volumes/ProjectData/TestPanelModel_tomahawk_abaqus_runs/TST-ART-02-HEAT-TRANFSER_001_shortened_name.in");
//    
//    mesh.prepare_for_use();
//    
//    // write the solution for visualization
//    libMesh::ExodusII_IO(mesh).write("mesh.exo");
//    
//    return 0;
//}
//
