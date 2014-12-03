
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




// C++ includes
#include <iostream>


// MAST includes


// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/abaqus_io.h"
#include "libmesh/exodusII_io.h"


int main(int argc, const char * argv[]) {
    
    // initialize the libMesh object
    libMesh::LibMeshInit              init(argc, argv);
    libMesh::ParallelMesh             mesh(init.comm());
    
    // initialize the mesh
    libMesh::AbaqusIO abaqus(mesh);
    abaqus.build_sidesets_from_nodesets = true;
    abaqus.read("/Volumes/ProjectData/TestPanelModel_tomahawk_abaqus_runs/TST-ART-02-HEAT-TRANFSER_001_shortened_name.in");
    
    mesh.prepare_for_use();
    
    // write the solution for visualization
    libMesh::ExodusII_IO(mesh).write("mesh.exo");
    
    return 0;
}

