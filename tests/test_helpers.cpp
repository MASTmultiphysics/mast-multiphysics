/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2020  Manav Bhatia and MAST authors
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

// We need access to the protected thermal_residual method to test it
// NOTE: Be careful with this, it could cause unexpected problems
#define protected public

#include "test_helpers.h"

#include "libmesh/point.h"
#include "libmesh/face_quad4.h"


#define pi 3.14159265358979323846


std::vector<double> TEST::eigen_matrix_to_std_vector(RealMatrixX M)
{
    std::vector<double> vec(M.data(), M.data()+M.rows()*M.cols());
    return vec;
}

Real TEST::get_shoelace_area(RealMatrixX X)
{
    Real true_volume = 0.0;
    uint n_nodes = X.cols();
    for (uint i=0; i<n_nodes-1; i++)
    {
        true_volume += X(0,i)*X(1,i+1) - X(0,i+1)*X(1,i);
    }
    true_volume += X(0,n_nodes-1)*X(1,0) - X(0,0)*X(1,n_nodes-1);
    true_volume = std::abs(true_volume);
    true_volume *= 0.5;
    return true_volume;
}

void TEST::approximate_internal_jacobian_with_finite_difference(
                                    MAST::StructuralElementBase& elem, 
                                    const RealVectorX& initial_elem_solution, 
                                    RealMatrixX& jacobian)
{
    Real delta = 3.4526698e-04; // sqrt(eps('single')) in MATLAB
    //Real delta = 1.1920929e-07; // eps('single') in MATLAB
    //Real delta = 1.490116119384766e-08; // sqrt(eps('double')) in MATLAB
    
    int n = jacobian.cols();
    RealMatrixX dummy = RealMatrixX::Zero(n,n);
    RealVectorX elem_solution = initial_elem_solution;    
    
    // 6th order central finite difference
    for (uint i=0; i<elem_solution.size(); i++)
    {
        elem_solution(i) = initial_elem_solution(i) + delta;
        elem.set_solution(elem_solution);
        RealVectorX residual_h = RealVectorX::Zero(n);
        elem.internal_residual(false, residual_h, dummy);
        
        elem_solution(i) = initial_elem_solution(i) + 2.0*delta;
        elem.set_solution(elem_solution);
        RealVectorX residual_2h = RealVectorX::Zero(n);
        elem.internal_residual(false, residual_2h, dummy);
        
        elem_solution(i) = initial_elem_solution(i) + 3.0*delta;
        elem.set_solution(elem_solution);
        RealVectorX residual_3h = RealVectorX::Zero(n);
        elem.internal_residual(false, residual_3h, dummy);
        
        elem_solution(i) = initial_elem_solution(i) - delta;
        elem.set_solution(elem_solution);
        RealVectorX residual_n = RealVectorX::Zero(n);
        elem.internal_residual(false, residual_n, dummy);
        
        elem_solution(i) = initial_elem_solution(i) - 2.0*delta;
        elem.set_solution(elem_solution);
        RealVectorX residual_2n = RealVectorX::Zero(n);
        elem.internal_residual(false, residual_2n, dummy);
        
        elem_solution(i) = initial_elem_solution(i) - 3.0*delta;
        elem.set_solution(elem_solution);
        RealVectorX residual_3n = RealVectorX::Zero(n);
        elem.internal_residual(false, residual_3n, dummy);
        
        elem_solution(i) = initial_elem_solution(i);
        
        //jacobian.col(i) = (residual_2n - 8.0*residual_n + 8.0*residual_h - residual_2h)/(12.0 * delta);
        
        jacobian.col(i) = (-residual_3n + 9.0*residual_2n - 45.0*residual_n + 45.0*residual_h - 9.0*residual_2h + residual_3h)/(60.0 * delta);
    }
}


void TEST::approximate_side_external_jacobian_with_finite_difference(
                                    MAST::StructuralElementBase& elem, 
                                    MAST::PhysicsDisciplineBase& discipline,
                                    const RealVectorX& initial_elem_solution, 
                                    RealMatrixX& jacobian)
{
    Real delta = 3.4526698e-04; // sqrt(eps('single')) in MATLAB
    
    int n = jacobian.cols();
    RealMatrixX dummy = RealMatrixX::Zero(n,n);
    RealVectorX elem_solution = initial_elem_solution;    
    
    // 4th order central finite difference
    for (uint i=0; i<elem_solution.size(); i++)
    {
        elem_solution(i) = initial_elem_solution(i) + delta;
        elem.set_solution(elem_solution);
        RealVectorX residual_h = RealVectorX::Zero(n);
        elem.side_external_residual(false, residual_h, dummy, dummy, discipline.side_loads());
        
        elem_solution(i) = initial_elem_solution(i) + 2.0*delta;
        elem.set_solution(elem_solution);
        RealVectorX residual_2h = RealVectorX::Zero(n);
        elem.side_external_residual(false, residual_2h, dummy, dummy, discipline.side_loads());
        
        elem_solution(i) = initial_elem_solution(i) - delta;
        elem.set_solution(elem_solution);
        RealVectorX residual_n = RealVectorX::Zero(n);
        elem.side_external_residual(false, residual_n, dummy, dummy, discipline.side_loads());
        
        elem_solution(i) = initial_elem_solution(i) - 2.0*delta;
        elem.set_solution(elem_solution);
        RealVectorX residual_2n = RealVectorX::Zero(n);
        elem.side_external_residual(false, residual_2n, dummy, dummy, discipline.side_loads());
        
        elem_solution(i) = initial_elem_solution(i);
        
        jacobian.col(i) = (residual_2n - 8.0*residual_n + 8.0*residual_h - residual_2h)/(12.0 * delta);
    }
}


void TEST::approximate_volume_external_jacobian_with_finite_difference(
                                    MAST::StructuralElementBase& elem, 
                                    MAST::PhysicsDisciplineBase& discipline,
                                    const RealVectorX& initial_elem_solution, 
                                    RealMatrixX& jacobian)
{
    Real delta = 3.4526698e-04; // sqrt(eps('single')) in MATLAB
    
    int n = jacobian.cols();
    RealMatrixX dummy = RealMatrixX::Zero(n,n);
    RealVectorX elem_solution = initial_elem_solution;    
    
    // 4th order central finite difference
    for (uint i=0; i<elem_solution.size(); i++)
    {
        elem_solution(i) = initial_elem_solution(i) + delta;
        elem.set_solution(elem_solution);
        RealVectorX residual_h = RealVectorX::Zero(n);
        elem.volume_external_residual(false, residual_h, dummy, dummy, discipline.volume_loads());
        
        elem_solution(i) = initial_elem_solution(i) + 2.0*delta;
        elem.set_solution(elem_solution);
        RealVectorX residual_2h = RealVectorX::Zero(n);
        elem.volume_external_residual(false, residual_2h, dummy, dummy, discipline.volume_loads());
        
        elem_solution(i) = initial_elem_solution(i) - delta;
        elem.set_solution(elem_solution);
        RealVectorX residual_n = RealVectorX::Zero(n);
        elem.volume_external_residual(false, residual_n, dummy, dummy, discipline.volume_loads());
        
        elem_solution(i) = initial_elem_solution(i) - 2.0*delta;
        elem.set_solution(elem_solution);
        RealVectorX residual_2n = RealVectorX::Zero(n);
        elem.volume_external_residual(false, residual_2n, dummy, dummy, discipline.volume_loads());
        
        elem_solution(i) = initial_elem_solution(i);
        
        jacobian.col(i) = (residual_2n - 8.0*residual_n + 8.0*residual_h - residual_2h)/(12.0 * delta);
    }
}


void TEST::approximate_inertial_jacobian_with_finite_difference(
                                    MAST::StructuralElementBase& elem, 
                                    const RealVectorX& initial_elem_accel, 
                                    RealMatrixX& jacobian)
{
    Real delta = 3.4526698e-04; // sqrt(eps('single')) in MATLAB
    //Real delta = 1.1920929e-07; // eps('single') in MATLAB
    //Real delta = 1.490116119384766e-08; // sqrt(eps('double')) in MATLAB
    
    int n = jacobian.cols();
    RealMatrixX dummy = RealMatrixX::Zero(n,n);
    RealVectorX elem_accel = initial_elem_accel;    
    
    // 6th order central finite difference
    for (uint i=0; i<elem_accel.size(); i++)
    {
        elem_accel(i) = initial_elem_accel(i) + delta;
        elem.set_acceleration(elem_accel);
        RealVectorX residual_h = RealVectorX::Zero(n);
        elem.inertial_residual(false, residual_h, dummy, dummy, dummy);
        
        elem_accel(i) = initial_elem_accel(i) + 2.0*delta;
        elem.set_acceleration(elem_accel);
        RealVectorX residual_2h = RealVectorX::Zero(n);
        elem.inertial_residual(false, residual_2h, dummy, dummy, dummy);
        
        elem_accel(i) = initial_elem_accel(i) + 3.0*delta;
        elem.set_acceleration(elem_accel);
        RealVectorX residual_3h = RealVectorX::Zero(n);
        elem.inertial_residual(false, residual_3h, dummy, dummy, dummy);
        
        elem_accel(i) = initial_elem_accel(i) - delta;
        elem.set_acceleration(elem_accel);
        RealVectorX residual_n = RealVectorX::Zero(n);
        elem.inertial_residual(false, residual_n, dummy, dummy, dummy);
        
        elem_accel(i) = initial_elem_accel(i) - 2.0*delta;
        elem.set_acceleration(elem_accel);
        RealVectorX residual_2n = RealVectorX::Zero(n);
        elem.inertial_residual(false, residual_2n, dummy, dummy, dummy);
        
        elem_accel(i) = initial_elem_accel(i) - 3.0*delta;
        elem.set_acceleration(elem_accel);
        RealVectorX residual_3n = RealVectorX::Zero(n);
        elem.inertial_residual(false, residual_3n, dummy, dummy, dummy);
        
        elem_accel(i) = initial_elem_accel(i);
        
        //jacobian.col(i) = (residual_2n - 8.0*residual_n + 8.0*residual_h - residual_2h)/(12.0 * delta);
        
        jacobian.col(i) = (-residual_3n + 9.0*residual_2n - 45.0*residual_n + 45.0*residual_h - 9.0*residual_2h + residual_3h)/(60.0 * delta);
    }
}


void TEST::approximate_thermal_jacobian_with_finite_difference(
                                    MAST::StructuralElementBase& elem, 
                                    const RealVectorX& initial_elem_solution, 
                                    RealMatrixX& jacobian,
                                    MAST::BoundaryConditionBase& thermal_bc)
{
    Real delta = 3.4526698e-04; // sqrt(eps('single')) in MATLAB
    //Real delta = 1.1920929e-07; // eps('single') in MATLAB
    //Real delta = 1.490116119384766e-08; // sqrt(eps('double')) in MATLAB
    
    int n = jacobian.cols();
    RealMatrixX dummy = RealMatrixX::Zero(n,n);
    RealVectorX elem_solution = initial_elem_solution;    
    
    // 6th order central finite difference
    for (uint i=0; i<elem_solution.size(); i++)
    {
        elem_solution(i) = initial_elem_solution(i) + delta;
        elem.set_solution(elem_solution);
        RealVectorX residual_h = RealVectorX::Zero(n);
        elem.thermal_residual(false, residual_h, dummy, thermal_bc);
        
        elem_solution(i) = initial_elem_solution(i) + 2.0*delta;
        elem.set_solution(elem_solution);
        RealVectorX residual_2h = RealVectorX::Zero(n);
        elem.thermal_residual(false, residual_2h, dummy, thermal_bc);
        
        elem_solution(i) = initial_elem_solution(i) + 3.0*delta;
        elem.set_solution(elem_solution);
        RealVectorX residual_3h = RealVectorX::Zero(n);
        elem.thermal_residual(false, residual_3h, dummy, thermal_bc);
        
        elem_solution(i) = initial_elem_solution(i) - delta;
        elem.set_solution(elem_solution);
        RealVectorX residual_n = RealVectorX::Zero(n);
        elem.thermal_residual(false, residual_n, dummy, thermal_bc);
        
        elem_solution(i) = initial_elem_solution(i) - 2.0*delta;
        elem.set_solution(elem_solution);
        RealVectorX residual_2n = RealVectorX::Zero(n);
        elem.thermal_residual(false, residual_2n, dummy, thermal_bc);
        
        elem_solution(i) = initial_elem_solution(i) - 3.0*delta;
        elem.set_solution(elem_solution);
        RealVectorX residual_3n = RealVectorX::Zero(n);
        elem.thermal_residual(false, residual_3n, dummy, thermal_bc);
        
        elem_solution(i) = initial_elem_solution(i);
        
        //jacobian.col(i) = (residual_2n - 8.0*residual_n + 8.0*residual_h - residual_2h)/(12.0 * delta);
        
        jacobian.col(i) = (-residual_3n + 9.0*residual_2n - 45.0*residual_n + 45.0*residual_h - 9.0*residual_2h + residual_3h)/(60.0 * delta);
    }
}


void TEST::transform_element(libMesh::MeshBase& mesh, const RealMatrixX X0,
                        Real shift_x, Real shift_y, Real shift_z,
                        Real scale_x, Real scale_y,
                        Real rotation_x, Real rotation_y, Real rotation_z,
                        Real shear_x, Real shear_y)
{    
    const uint n_nodes = X0.cols();
    
    RealMatrixX Vi = RealMatrixX::Ones(1,n_nodes);
    
    //std::cout << "X0:\n" << X0 << "\n" << std::endl;
    
    // Vectors from element centroid to its nodes
    RealMatrixX Vc = RealMatrixX::Zero(3,n_nodes);
    Vc.row(0) = X0.row(0) - X0.row(0).mean()*Vi;
    Vc.row(1) = X0.row(1) - X0.row(1).mean()*Vi;
    Vc.row(2) = X0.row(2) - X0.row(2).mean()*Vi;
    
    //std::cout << "Vc:\n" << Vc << "\n" << std::endl;
    
    // Scaling Matrix
    RealMatrixX S = RealMatrixX::Identity(3,3);
    S(0,0) = scale_x;
    S(1,1) = scale_y;
    
    //std::cout << "S:\n" << S << "\n" << std::endl;
    
    // Scale the element's vectors
    RealMatrixX V = S * Vc;
    
    //std::cout << "V:\n" << V << "\n" << std::endl;
    
    
    // Shear Matrix
    RealMatrixX K = RealMatrixX::Identity(3,3);
    K(0,1) = shear_x;
    K(1,0) = shear_y;
    
    // Shear the element's vectors
    V = K * V;
    
    // Rotation Matrices
    RealMatrix3 Rz;
    Real theta = rotation_z*pi/180.0;
    Rz <<  cos(theta),    -sin(theta),         0.0,
           sin(theta),     cos(theta),         0.0,
             0.0,             0.0,             1.0;
    //std::cout << "Rz:\n" << Rz << "\n" << std::endl;
             
    RealMatrix3 Ry;
    theta = rotation_y*pi/180.0;
    Ry <<  cos(theta),        0.0,         -sin(theta),
             0.0,             1.0,             0.0,
           sin(theta),        0.0,          cos(theta);
    //std::cout << "Ry:\n" << Ry << "\n" << std::endl;
           
    RealMatrix3 Rx;
    theta = rotation_x*pi/180.0;
    Rx <<    1.0,             0.0,                 0.0,
             0.0,          cos(theta),         -sin(theta),
             0.0,          sin(theta),          cos(theta);
    //std::cout << "Rx:\n" << Rx << "\n" << std::endl;
             
    // Rotate the element's vectors about the centroid
    V = Rz * V;
    V = Ry * V;
    V = Rx * V;
    
    //std::cout << "V:\n" << V << "\n" << std::endl;
    
    RealMatrixX X = RealMatrixX::Zero(3,n_nodes);
    
    //std::cout << "xc:\n" << X0.row(0).mean()*Vi << "\n" << std::endl;
    //std::cout << "yc:\n" << X0.row(1).mean()*Vi << "\n" << std::endl;
    //std::cout << "zc:\n" << X0.row(2).mean()*Vi << "\n" << std::endl;
    
    X.row(0) = X0.row(0).mean()*Vi + V.row(0);
    X.row(1) = X0.row(1).mean()*Vi + V.row(1);
    X.row(2) = X0.row(2).mean()*Vi + V.row(2);
    
    // Shift the element
    X.row(0) += shift_x*RealMatrixX::Ones(1,n_nodes);
    X.row(1) += shift_y*RealMatrixX::Ones(1,n_nodes);
    X.row(2) += shift_z*RealMatrixX::Ones(1,n_nodes);
    
    //std::cout << "Transformed X:\n" << X << "\n" << std::endl;
    
    
    // Update the element with the new node Coordinates
    for (int i=0; i<n_nodes; i++)
    {
        (*mesh.node_ptr(i)) = libMesh::Point(X(0,i), X(1,i), X(2,i));
    }
}
