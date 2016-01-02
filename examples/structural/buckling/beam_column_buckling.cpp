
// C++ includes
#include <iostream>


// MAST includes
#include "examples/structural/buckling/beam_column_buckling.h"
#include "base/nonlinear_system.h"
#include "elasticity/structural_system_initialization.h"
#include "elasticity/structural_discipline.h"
#include "elasticity/stress_output_base.h"
#include "property_cards/solid_1d_section_element_property_card.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "property_cards/isotropic_material_property_card.h"
#include "elasticity/structural_element_base.h"
#include "elasticity/structural_modal_eigenproblem_assembly.h"
#include "boundary_condition/dirichlet_boundary_condition.h"
#include "driver/driver_base.h"


// libMesh includes
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/numeric_vector.h"

extern libMesh::LibMeshInit* __init;


MAST::BeamColumnBucklingAnalysis::BeamColumnBucklingAnalysis() {
    
    
    // create the mesh
    _mesh       = new libMesh::SerialMesh(__init->comm());
    
    // initialize the mesh with one element
    libMesh::MeshTools::Generation::build_line(*_mesh, 10, 0, 10);
    _mesh->prepare_for_use();
    
    // create the equation system
    _eq_sys    = new  libMesh::EquationSystems(*_mesh);
    
    // create the libmesh system
    _sys       = &(_eq_sys->add_system<MAST::NonlinearSystem>("structural"));
    _sys->set_eigenproblem_type(libMesh::GHEP);

    // FEType to initialize the system
    libMesh::FEType fetype (libMesh::FIRST, libMesh::LAGRANGE);
    
    // initialize the system to the right set of variables
    _structural_sys = new MAST::StructuralSystemInitialization(*_sys,
                                                               _sys->name(),
                                                               fetype);
    _discipline     = new MAST::StructuralDiscipline(*_eq_sys);
    
    
    // create and add the boundary condition and loads
    _dirichlet_left = new MAST::DirichletBoundaryCondition;
    _dirichlet_right= new MAST::DirichletBoundaryCondition;
    std::vector<unsigned int> constrained_vars(4);
    constrained_vars[0] = 0;  // u
    constrained_vars[1] = 1;  // v
    constrained_vars[2] = 2;  // w
    constrained_vars[3] = 3;  // tx
    _dirichlet_left->init (0, constrained_vars);
    _dirichlet_right->init(1, constrained_vars);
    _discipline->add_dirichlet_bc(0, *_dirichlet_left);
    _discipline->add_dirichlet_bc(1, *_dirichlet_right);
    _discipline->init_system_dirichlet_bc(*_sys);
    
    // initialize the equation system
    _eq_sys->init();
    
    _sys->eigen_solver->set_position_of_spectrum(libMesh::LARGEST_MAGNITUDE);

    // create the property functions and add them to the
    
    _thy             = new MAST::Parameter("thy", 0.06);
    _thz             = new MAST::Parameter("thz", 0.02);
    _rho             = new MAST::Parameter("rho",2.8e3);
    _E               = new MAST::Parameter("E",  72.e9);
    _nu              = new MAST::Parameter("nu",  0.33);
    _zero            = new MAST::Parameter("zero",  0.);
    
    
    
    // prepare the vector of parameters with respect to which the sensitivity
    // needs to be benchmarked
    _params_for_sensitivity.push_back(_E);
    _params_for_sensitivity.push_back(_nu);
    _params_for_sensitivity.push_back(_thy);
    _params_for_sensitivity.push_back(_thz);
    
    
    
    _thy_f           = new MAST::ConstantFieldFunction("hy",     *_thy);
    _thz_f           = new MAST::ConstantFieldFunction("hz",     *_thz);
    _rho_f           = new MAST::ConstantFieldFunction("rho",    *_rho);
    _E_f             = new MAST::ConstantFieldFunction("E",      *_E);
    _nu_f            = new MAST::ConstantFieldFunction("nu",     *_nu);
    _hyoff_f         = new MAST::ConstantFieldFunction("hy_off", *_zero);
    _hzoff_f         = new MAST::ConstantFieldFunction("hz_off", *_zero);
    
    // create the material property card
    _m_card         = new MAST::IsotropicMaterialPropertyCard;
    
    // add the material properties to the card
    _m_card->add(*_rho_f);
    _m_card->add(*_E_f);
    _m_card->add(*_nu_f);
    
    // create the element property card
    _p_card         = new MAST::Solid1DSectionElementPropertyCard;
    
    // tell the card about the orientation
    libMesh::Point orientation;
    orientation(1) = 1.;
    _p_card->y_vector() = orientation;
    
    // add the section properties to the card
    _p_card->add(*_thy_f);
    _p_card->add(*_thz_f);
    _p_card->add(*_hyoff_f);
    _p_card->add(*_hzoff_f);
    
    // tell the section property about the material property
    _p_card->set_material(*_m_card);
    
    _p_card->init();
    
    _discipline->set_property_for_subdomain(0, *_p_card);
}







MAST::BeamColumnBucklingAnalysis::~BeamColumnBucklingAnalysis() {
    
    delete _m_card;
    delete _p_card;
    
    delete _dirichlet_left;
    delete _dirichlet_right;
    
    delete _thy_f;
    delete _thz_f;
    delete _rho_f;
    delete _E_f;
    delete _nu_f;
    delete _hyoff_f;
    delete _hzoff_f;
    
    delete _thy;
    delete _thz;
    delete _rho;
    delete _E;
    delete _nu;
    delete _zero;
    
    
    
    
    delete _eq_sys;
    delete _mesh;
    
    delete _discipline;
    delete _structural_sys;
}




const libMesh::NumericVector<Real>&
MAST::BeamColumnBucklingAnalysis::solve() {
    
    
    // create the nonlinear assembly object
    MAST::StructuralModalEigenproblemAssembly   assembly;
    assembly.set_exchange_A_and_B_matrices(true);
    _sys->initialize_condensed_dofs(*_discipline);

    const unsigned int n_eig = 4;
    _eq_sys->parameters.set<unsigned int>("eigenpairs")    = n_eig;
    _eq_sys->parameters.set<unsigned int>("basis vectors") = n_eig*3;
    
    
    assembly.attach_discipline_and_system(*_discipline, *_structural_sys);
    _sys->solve();
    assembly.clear_discipline_and_system();
    
    // Get the number of converged eigen pairs.
    unsigned int nconv = std::min(_sys->get_n_converged(), n_eig);
    
    for (unsigned int i=0; i<nconv; i++)
    {
        std::ostringstream file_name;
        
        // We write the file in the ExodusII format.
        file_name << "out_"
        << std::setw(3)
        << std::setfill('0')
        << std::right
        << i
        << ".exo";
        
        // now write the eigenvlaues
        std::pair<Real, Real> val = _sys->get_eigenpair(i);
        std::complex<Real> eigval = std::complex<Real>(val.first, val.second);
        eigval = 1./eigval;
        libMesh::out
        << std::setw(35) << std::fixed << std::setprecision(15)
        << eigval.real() << std::endl;
        
        // We write the file in the ExodusII format.
        libMesh::ExodusII_IO(*_mesh).write_equation_systems(file_name.str(),
                                                            *_eq_sys);
    }

    
    return *(_sys->solution);
}





const libMesh::NumericVector<Real>&
MAST::BeamColumnBucklingAnalysis::sensitivity_solve(MAST::Parameter& p) {
    
    
//    // create the nonlinear assembly object
//    MAST::StructuralNonlinearAssembly   assembly;
//    
//    // now solve the system
//    MAST::Driver::sensitivity_solution(*_discipline,
//                                       *_structural_sys,
//                                       assembly,
//                                       p);
//    
//    // write the solution for visualization
//    //libMesh::ExodusII_IO(mesh).write_equation_systems("mesh.exo", eq_sys);
//    
    return _sys->get_sensitivity_solution(0);
}



