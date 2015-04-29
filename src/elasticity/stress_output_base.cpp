

// MAST includes
#include "elasticity/stress_output_base.h"



MAST::StressStrainOutputBase::StressStrainOutputBase():
MAST::OutputFunctionBase(),
_if_evaluate_stress(true),
_if_evaluate_strain(false) {
    
}



MAST::StressStrainOutputBase::~StressStrainOutputBase() {
    
}



void
MAST::StressStrainOutputBase::set_evaluate_stress(bool val) {
    _if_evaluate_stress = val;
}



void
MAST::StressStrainOutputBase::set_evaluate_strain(bool val) {
    _if_evaluate_strain = val;
}


bool
MAST::StressStrainOutputBase::if_evaluate_stress() const {
    return _if_evaluate_stress;
}



bool
MAST::StressStrainOutputBase::if_evaluate_strain() const {
    return _if_evaluate_strain;
}



void
MAST::StressStrainOutputBase::
add_stress_at_qp_location(const libMesh::Point& quadrature_pt,
                          const libMesh::Point& physical_pt,
                          const RealVectorX& vec) {
    
}



void
MAST::StressStrainOutputBase::
add_strain_at_qp_location(const libMesh::Point& quadrature_pt,
                          const libMesh::Point& physical_pt,
                          const RealVectorX& vec) {
    
}

