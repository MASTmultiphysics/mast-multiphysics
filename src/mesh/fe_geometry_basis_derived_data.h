
#ifndef _mast_fe_geometry_basis_derived_data_h_
#define _mast_fe_geometry_basis_derived_data_h_

// MAST includes
#include "mesh/fe_shape_data_base.h"


namespace MAST {

/*! This provides the derivative of shape functions when the FE basis also forms the basis for geometry interpolation used to interpolate
 nodal locations. Typically, Lagrange shape functoins are used for this purpose */
template <typename BasisScalarType, typename NodalScalarType, typename Traits, typename ContextType>
class FEGeometryBasisDerivedData:
public MAST::FEShapeDataBase<BasisScalarType, NodalScalarType, Traits, ContextType> {
    
public:
    
    FEGeometryBasisDerivedData(const std::string& nm):
    MAST::FEShapeDataBase<BasisScalarType, NodalScalarType, Traits, ContextType>(nm) {}
    virtual ~FEGeometryBasisDerivedData() {}
    
protected:

};



/*! This provides the derivative of shape functions when the FE basis also forms the basis for geometry interpolation used to interpolate
 nodal locations. Typically, Lagrange shape functoins are used for this purpose */
template <typename BasisScalarType, typename NodalScalarType, typename ContextType>
class FEGeometryBasisDerivedData<BasisScalarType, NodalScalarType, EigenTraits, ContextType>:
public MAST::FEShapeDataBase<BasisScalarType, NodalScalarType, EigenTraits, ContextType> {
    
public:
            
    FEGeometryBasisDerivedData(const std::string& nm):
    MAST::FEShapeDataBase<BasisScalarType, NodalScalarType, EigenTraits, ContextType>(nm) {}
    
    virtual ~FEGeometryBasisDerivedData() {}
    
    virtual inline void reinit(const ContextType& c) override {
        
        uint_type
        nq  = this->n_q_points();
        const int_type
        d   = this->spatial_dim();
        const uint_type
        n_nodes = c.n_nodes();
        
        // for this class the number of basis functions should be equal to the number
        // of nodes
        libmesh_assert_equal_to(n_nodes, this->_fe_basis->n_basis());
        
        if (this->_compute_xyz) {
            
            this->_node_coord  = EigenMatrix<NodalScalarType>::type::Zero(n_nodes, d);
            this->_xyz         = EigenMatrix<NodalScalarType>::type::Zero(nq, d);

            // get the nodal locations
            for (uint_type i=0; i<n_nodes; i++)
                for (uint_type j=0; j<d; j++)
                    this->_node_coord(i, j) = c.nodal_coord(i, j);

            // quadrature point coordinates
            for (uint_type i=0; i<nq; i++)
                for (uint_type k=0; k<n_nodes; k++)
                    for (uint_type j=0; j<d; j++)
                        this->_xyz(i, j) += this->_fe_basis->phi(i, k) * this->_node_coord(k, j);
        }
        
        
        if (this->_compute_Jac) {
            
            this->_dx_dxi      = EigenMatrix<NodalScalarType>::type::Zero(nq, d*d);

            // quadrature point spatial coordinate derivatives dx/dxi
            for (uint_type i=0; i<nq; i++) {
                
                Eigen::Map<typename EigenMatrix<NodalScalarType>::type>
                dxdxi(this->_dx_dxi.row(i).data(), d, d);

                for (uint_type l=0; l<n_nodes; l++)
                    for (uint_type j=0; j<d; j++)
                        for (uint_type k=0; k<d; k++)
                            dxdxi(j, k) += this->_fe_basis->dphi_dxi(i, l, k) * this->_node_coord(l, j);
            }
        }
        
        
        if (this->_compute_detJ) {

            this->_detJ        = EigenVector<NodalScalarType>::type::Zero(nq);
            
            for (uint_type i=0; i<nq; i++) {
                
                Eigen::Map<typename EigenMatrix<NodalScalarType>::type>
                dxdxi(this->_dx_dxi.row(i).data(), d, d);

                // determinant of dx_dxi
                this->_detJ(i) = dxdxi.determinant();
            }
        }

        
        if (this->_compute_JxW) {

            this->_detJxW      = EigenVector<NodalScalarType>::type::Zero(nq);
            
            for (uint_type i=0; i<nq; i++) {
                
                Eigen::Map<typename EigenMatrix<NodalScalarType>::type>
                dxdxi(this->_dx_dxi.row(i).data(), d, d);

                // determinant times weight
                this->_detJxW(i) = this->_detJ(i) * this->_fe_basis->qp_weight(i);
            }
        }

        
        if (this->_compute_Jac_inv) {
            
            this->_dxi_dx      = EigenMatrix<NodalScalarType>::type::Zero(nq, d*d);

            for (uint_type i=0; i<nq; i++) {

                // quadrature point spatial coordinate derivatives dx/dxi
                Eigen::Map<typename EigenMatrix<NodalScalarType>::type>
                dxdxi(this->_dx_dxi.row(i).data(), d, d),
                dxidx(this->_dxi_dx.row(i).data(), d, d);
                
                // compute dx/dxi
                dxidx = dxdxi.inverse();
            }
        }
        
        
        if (this->_compute_dphi_dx) {

            this->_dphi_dx     = EigenMatrix<NodalScalarType>::type::Zero(nq, d*d*this->n_basis());
            
            for (uint_type i=0; i<nq; i++) {
                
                // quadrature point spatial coordinate derivatives dx/dxi
                Eigen::Map<typename EigenMatrix<NodalScalarType>::type>
                dxidx (this->_dxi_dx.row(i).data(), d, d),
                dphidx(this->_dphi_dx.row(i).data(), d, n_nodes);
                
                for (uint_type l=0; l<n_nodes; l++)
                    for (uint_type j=0; j<d; j++)
                        for (uint_type k=0; k<d; k++)
                            dphidx(j, l) += this->_fe_basis->dphi_dxi(i, l, k) * dxidx(k, j);
            }
        }
    }

    virtual inline void reinit_for_side(const ContextType& c, uint_type s) override { }
    
    
protected:

};

}

#endif // _mast_fe_geometry_basis_derived_data_h_
