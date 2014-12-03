
// C++ includes
#include <iostream>


// MAST includes
#include "heat_conduction/heat_conduction_discipline.h"
#include "heat_conduction/heat_conduction_system_initialization.h"
#include "heat_conduction/heat_conduction_nonlinear_assembly.h"
#include "base/physics_discipline_base.h"
#include "base/output_assembly_base.h"
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
#include "libmesh/mesh_function.h"
#include "libmesh/mesh_base.h"
#include "libmesh/node.h"


class TemperatureOutput:
public MAST::OutputAssemblyBase {
    
public:
    
    TemperatureOutput(MAST::SystemInitialization& system,
                      const libMesh::Point p):
    MAST::OutputAssemblyBase(),
    _p(p),
    _sys(system.system()),
    _var_num(system.vars()[0]) {
        
        _sys.qoi.resize(1);
        
        libMesh::MeshFunction *rval =
        new libMesh::MeshFunction(_sys.get_equation_systems(),
                                  *_sys.solution,
                                  _sys.get_dof_map(),
                                  system.vars());
        rval->init();
        
        _func.reset(rval);
    }
    
    virtual ~TemperatureOutput() {}
    
    
    virtual void qoi (const libMesh::QoISet& qoi_indices) {
        // temperature at the center of the block
        DenseRealVector v(1);
        (*_func)(_p, 0., v);
        _sys.qoi[0] = v(0);
    }
    
    
    virtual void qoi_derivative (const libMesh::QoISet& qoi_indices,
                                 bool include_liftfunc,
                                 bool apply_constraints) {
        
        libMesh::NumericVector<Real>& deriv = _sys.get_adjoint_rhs();
        deriv.zero();
        
        // if there is a node at the given point, then add a unit value
        // for its dof.
        libMesh::MeshBase& m = _sys.get_mesh();
        libMesh::MeshBase::const_node_iterator
        it    = m.local_nodes_begin(),
        end   = m.local_nodes_end();
        
        for ( ; it != end; it++)
            if (_p == **it)
                deriv.set((*it)->dof_number(_sys.number(), _var_num, 0), 1.);
        
        deriv.close();
    }
    
    
    virtual bool
    qoi_parameter_sensitivity (const libMesh::QoISet& qoi_indices,
                               const libMesh::ParameterVector& parameters,
                               const unsigned int p,
                               std::vector<Real>& partialq_partialp) {
        
        // assumed zero for this parameter
        partialq_partialp[0] = 0.;
        
        return true;
    }

protected:
    
    virtual std::auto_ptr<MAST::ElementBase>
    _build_elem(const libMesh::Elem& elem) {
        libmesh_not_implemented();
    }

    libMesh::Point _p;
    libMesh::System& _sys;
    const unsigned int _var_num;
    std::auto_ptr<libMesh::MeshFunction> _func;
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
    MAST::Parameter k  (  "k", 2.0);   // thermal conductance
    MAST::Parameter c  (  "c", 2.0);   // thermal capacitance
    MAST::Parameter rho("rho", 2.0); // material density
    
    
    // define material property based on the parameter
    MAST::ConstantFieldFunction  k_th("k_th", k);
    MAST::ConstantFieldFunction  c_th("c_th", c);
    
    
    // add them to a material property card
    MAST::IsotropicMaterialPropertyCard       materials(1);
    materials.add(k_th);
    materials.add(c_th);
    
    
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
    libMesh::ExodusII_IO(mesh).write_equation_systems("primal.exo", eq_sys);
    
    // Evaluate the output function
    TemperatureOutput t_out(heat_cond_sys,
                            libMesh::Point(width/2., height/2., 0.));
    
    // Evaluate the output
    MAST::Driver::output_evaluation(heat_cond, heat_cond_sys, t_out);
    std::cout << "Output Function: " << sys.qoi[0] << std::endl;

    // Calculate the adjoint solution
    MAST::Driver::adjoint_solution(heat_cond, heat_cond_sys, assembly, t_out);
    
    // write the adjoint solution for visualization
    sys.solution->swap(sys.get_adjoint_solution());
    libMesh::ExodusII_IO(mesh).write_equation_systems("adjoint.exo", eq_sys);
    
    return 0;
}

