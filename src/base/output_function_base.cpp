
// MAST includes
#include "base/output_function_base.h"

MAST::OutputFunctionBase::OutputFunctionBase() {
    
}



MAST::OutputFunctionBase::~OutputFunctionBase() {
    
}



void
MAST::OutputFunctionBase::
set_qp_for_evaluation(const std::vector<libMesh::Point>& pts) {
    
    // make sure that some points were specified
    libmesh_assert(pts.size());
    
    _eval_mode = MAST::SPECIFIED_QP;
    _eval_points = pts;
}



const std::vector<libMesh::Point>&
MAST::OutputFunctionBase::get_qp_for_evaluation() const {

    // make sure that the data was provide
    libmesh_assert(_eval_mode == MAST::SPECIFIED_QP);
    return _eval_points;
}


