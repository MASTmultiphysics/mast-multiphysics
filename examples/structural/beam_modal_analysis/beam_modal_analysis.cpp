
// C++ includes
#include <iostream>


// MAST includes
#include "examples/structural/beam_modal_analysis/beam_modal_analysis.h"
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


MAST::BeamModalAnalysis::BeamModalAnalysis() {
    
    
    // create the mesh
    _mesh       = new libMesh::SerialMesh(__init->comm());
    _length     = 10.;
    
    // initialize the mesh with one element
    libMesh::MeshTools::Generation::build_line(*_mesh, 50, 0, _length);
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
    _sys->set_exchange_A_and_B(true);
    _sys->set_n_requested_eigenvalues(3);
    
    // create the property functions and add them to the
    
    _thy             = new MAST::Parameter("thy", 0.06);
    _thz             = new MAST::Parameter("thz", 1.00);
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







MAST::BeamModalAnalysis::~BeamModalAnalysis() {
    
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




MAST::Parameter*
MAST::BeamModalAnalysis::get_parameter(const std::string &nm) {
    
    MAST::Parameter *rval = NULL;
    
    // look through the vector of parameters to see if the name is available
    std::vector<MAST::Parameter*>::iterator
    it   =  _params_for_sensitivity.begin(),
    end  =  _params_for_sensitivity.end();
    
    bool
    found = false;
    
    for ( ; it != end; it++) {
        
        if (nm == (*it)->name()) {
            rval    = *it;
            found   = true;
        }
    }
    
    // if the param was not found, then print the message
    if (!found) {
        std::cout
        << std::endl
        << "Parameter not found by name: " << nm << std::endl
        << "Valid names are: "
        << std::endl;
        for (it = _params_for_sensitivity.begin(); it != end; it++)
            std::cout << "   " << (*it)->name() << std::endl;
        std::cout << std::endl;
    }
    
    return rval;
}



void
MAST::BeamModalAnalysis::solve(std::vector<Real>& eig,
                               bool if_write_output) {
    
    
    // create the nonlinear assembly object
    MAST::StructuralModalEigenproblemAssembly   assembly;
    _sys->initialize_condensed_dofs(*_discipline);
    
    assembly.attach_discipline_and_system(*_discipline, *_structural_sys);
    _sys->eigenproblem_solve();
    assembly.clear_discipline_and_system();
    
    // Get the number of converged eigen pairs.
    unsigned int
    nconv = std::min(_sys->get_n_converged_eigenvalues(),
                     _sys->get_n_requested_eigenvalues());
    eig.resize(nconv);
    
    for (unsigned int i=0; i<nconv; i++) {
        
        std::ostringstream file_name;
        
        // We write the file in the ExodusII format.
        file_name << "out_"
        << std::setw(3)
        << std::setfill('0')
        << std::right
        << i
        << ".exo";
        
        // now write the eigenvlau
        Real
        re = 0.,
        im = 0.;
        _sys->get_eigenpair(i, re, im, *_sys->solution);
        eig[i] = re;
        
        libMesh::out
        << std::setw(35) << std::fixed << std::setprecision(15)
        << re << std::endl;
        
        if (if_write_output) {
            
            std::cout
            << "Writing mode " << i << " to : "
            << file_name.str() << std::endl;
            

            // We write the file in the ExodusII format.
            libMesh::ExodusII_IO(*_mesh).write_equation_systems(file_name.str(),
                                                                *_eq_sys);
        }
    }
}





void
MAST::BeamModalAnalysis::sensitivity_solve(MAST::Parameter& p,
                                           std::vector<Real>& eig) {

    _discipline->add_parameter(p);

    // Get the number of converged eigen pairs.
    unsigned int
    nconv = std::min(_sys->get_n_converged_eigenvalues(),
                     _sys->get_n_requested_eigenvalues());
    eig.resize(nconv);

    libMesh::ParameterVector params;
    params.resize(1);
    params[0]  =  p.ptr();

    // create the nonlinear assembly object
    MAST::StructuralModalEigenproblemAssembly   assembly;
    assembly.attach_discipline_and_system(*_discipline, *_structural_sys);
    _sys->eigenproblem_sensitivity_solve(params, eig);
    assembly.clear_discipline_and_system();
    
    _discipline->remove_parameter(p);
}


