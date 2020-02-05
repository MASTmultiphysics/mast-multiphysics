/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2019  Manav Bhatia
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


//#define BOOST_TEST_DYN_LINK
//#define BOOST_TEST_MODULE InviscidFluxEigensystemDerivative
//#include <boost/test/unit_test.hpp>
//
//// Test includes
//#include "fluid/fluid_elem_initialization.h"
//
//
//
//BOOST_FIXTURE_TEST_SUITE(L1Derivative, BuildFluidElem)
//
//
//
//// sensitivity wrt rho is all zero. Hence, the first one tested is wrt u1
//BOOST_AUTO_TEST_CASE(L1U1DerivativeGeneralN)
//{
//    DenseRealMatrix eig(5,5), mat(5,5), mat_inv(5,5);
//    RealMatrixX eig0(5,5), mat0(5,5), mat0_inv(5,5);
//    
//    // initialize the analytical derivative for comparison
//    mat0 <<
//    121.74913841611915, 0.5933908290969268, 0.2373563316387707, 875.6885867287625, -655.6885867287625, -1, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
//    mat0_inv <<
//    0, 0, 0, 0, 0, -3.0212393123659325e-6, 0.00003323363243602526, 0.00006646726487205052, 1.5106196561829662e-6, 1.5106196561829662e-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.00033233632436025255, -0.4176517890908552, -0.16121159632760163, -0.0002974668013716487, 0.0006298031257319012;
//    
//    libMesh::Point n(sqrt(1.-.25-.04),0.5,0.2);
//
//    _fluid_elem->calculate_advection_left_eigenvector_and_inverse_u1_derivative_for_normal(*_primitive_sol,
//                                                                                           n,
//                                                                                           eig,
//                                                                                           mat,
//                                                                                           mat_inv);
//    //BOOST_CHECK(compare(eig0, eig));
//    BOOST_CHECK(compare(mat0, mat));
//    BOOST_CHECK(compare(mat0_inv, mat_inv));
//}
//
//
//// test wrt u2
//BOOST_AUTO_TEST_CASE(L1U2DerivativeGeneralN)
//{
//    DenseRealMatrix eig(5,5), mat(5,5), mat_inv(5,5);
//    RealMatrixX eig0(5,5), mat0(5,5), mat0_inv(5,5);
//    
//    // initialize the analytical derivative for comparison
//    mat0 <<
//    54.27299120066195, -1, 0, 465.35258530903457, -443.3525853090345, -0.5933908290969268, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
//    mat0_inv <<
//    0, 3.0212393123659325e-6, 0, 0, 0, 0, 0.00033233632436025255, 0, 0, 0, -3.0212393123659325e-6, 0.00006646726487205052, 0.00006646726487205052, 1.5106196561829662e-6, 1.5106196561829662e-6, 0, 0.00006646726487205052, 0, 0, 0, -0.00003323363243602526, 0.7695579926886009, -0.09926886008640748, -0.0002584999192022985, 0.0002917335516383238;
//    
//    libMesh::Point n(sqrt(1.-.25-.04),0.5,0.2);
//    _fluid_elem->calculate_advection_left_eigenvector_and_inverse_u2_derivative_for_normal(*_primitive_sol,
//                                                                                           n,
//                                                                                           eig,
//                                                                                           mat,
//                                                                                           mat_inv);
//    
//    //BOOST_CHECK(compare(eig0, eig));
//    BOOST_CHECK(compare(mat0, mat));
//    BOOST_CHECK(compare(mat0_inv, mat_inv));
//}
//
//
//// test wrt u3
//BOOST_AUTO_TEST_CASE(L1U3DerivativeGeneralN)
//{
//    DenseRealMatrix eig(5,5), mat(5,5), mat_inv(5,5);
//    RealMatrixX eig0(5,5), mat0(5,5), mat0_inv(5,5);
//    
//    // initialize the analytical derivative for comparison
//    mat0 <<
//    4.109196480264778, 0, -1, 203.74103412361382, -159.74103412361384, -0.2373563316387707, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0;
//    mat0_inv <<
//    0, 0, 3.0212393123659325e-6, 0, 0, 0, 0, 0.00033233632436025255, 0, 0, 0, 0, 0.00003323363243602526, 0, 0, -3.0212393123659325e-6, 0.00003323363243602526, 0.00013293452974410104, 1.5106196561829662e-6, 1.5106196561829662e-6, -0.00006646726487205052, -0.09926886008640748, 0.9806547025589896, -0.00007681306173209923, 0.00014328032660414973;
//    
//    libMesh::Point n(sqrt(1.-.25-.04),0.5,0.2);
//    _fluid_elem->calculate_advection_left_eigenvector_and_inverse_u3_derivative_for_normal(*_primitive_sol,
//                                                                                           n,
//                                                                                           eig,
//                                                                                           mat,
//                                                                                           mat_inv);
//    
//    //BOOST_CHECK(compare(eig0, eig));
//    BOOST_CHECK(compare(mat0, mat));
//    BOOST_CHECK(compare(mat0_inv, mat_inv));
//}
//
//
//
//// test wrt u3
//BOOST_AUTO_TEST_CASE(L1TDerivativeGeneralN)
//{
//    DenseRealMatrix eig(5,5), mat(5,5), mat_inv(5,5);
//    RealMatrixX eig0(5,5), mat0(5,5), mat0_inv(5,5);
//    
//    // initialize the analytical derivative for comparison
//    mat0 <<
//    -1003, 0, 0, 141.24534201406482, -141.24534201406482, 0, 0, 0, -1.1601342223163071, 1.1601342223163071, 0, 0, 0, -0.6884130080439917, 0.6884130080439917, 0, 0, 0, -0.2753652032175967, 0.2753652032175967, 0, 0, 0, 0, 0;
//    mat0_inv <<
//    9.15527064353313e-9, -1.0070797707886442e-7, -2.0141595415772884e-7, -4.577635321766565e-9, -4.577635321766565e-9, 1.0070797707886442e-6, -0.000011077877478675086, -0.000022155754957350172, 1.9893733210836735e-7, -1.2060171028970112e-6, 1.0070797707886442e-7, -1.1077877478675086e-6, -2.2155754957350172e-6, 3.6648954997619075e-7, -4.671975270550552e-7, 2.0141595415772884e-7, -2.2155754957350172e-6, -4.4311509914700345e-6, 6.602943832738479e-8, -2.6744539248511365e-7, 0.0000581588567630442, -0.0006397474243934863, -0.0012794948487869725, 0.00005644656760638308, -0.00011460542436942726;
//    
//    libMesh::Point n(sqrt(1.-.25-.04),0.5,0.2);
//    _fluid_elem->calculate_advection_left_eigenvector_and_inverse_T_derivative_for_normal(*_primitive_sol,
//                                                                                          n,
//                                                                                          eig,
//                                                                                          mat,
//                                                                                          mat_inv);
//    
//    //BOOST_CHECK(compare(eig0, eig));
//    BOOST_CHECK(compare(mat0, mat));
//    BOOST_CHECK(compare(mat0_inv, mat_inv));
//}
//
//
//BOOST_AUTO_TEST_SUITE_END()
//
//
//
//BOOST_FIXTURE_TEST_SUITE(L2Derivative, BuildFluidElem)
//
//
//// sensitivity wrt rho is all zero. Hence, the first one tested is wrt u1
//BOOST_AUTO_TEST_CASE(L2U1Derivative)
//{
//    DenseRealMatrix eig(5,5), mat(5,5), mat_inv(5,5);
//    RealMatrixX eig0(5,5), mat0(5,5), mat0_inv(5,5);
//    
//    // initialize the analytical derivative for comparison
//    mat0 <<
//    -1, -103.47270087993381, 0, 564.3525853090346, -344.3525853090345, 0, 0, 0, -1, -1, 0, -0.5933908290969268, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
//    mat0_inv <<
//    3.0212393123659325e-6, 0, 0, 0, 0, 0.0006646726487205051, -3.0212393123659325e-6, 0.00006646726487205052, 1.5106196561829662e-6, 1.5106196561829662e-6, 0.00003323363243602526, 0, 0, 0, 0, 0.00006646726487205052, 0, 0, 0, 0, 0.8057494184114324, -0.00033233632436025255, -0.09268860086407446, -0.00010894857324018491, 0.0004412848976004374;
//    
//    libMesh::Point n(0.5,sqrt(1.-.5*.5-.2*.2),0.2);
//    _fluid_elem->calculate_advection_left_eigenvector_and_inverse_u1_derivative_for_normal(*_primitive_sol,
//                                                                                           n,
//                                                                                           eig,
//                                                                                           mat,
//                                                                                           mat_inv);
//    
//    //BOOST_CHECK(compare(eig0, eig));
//    BOOST_CHECK(compare(mat0, mat));
//    BOOST_CHECK(compare(mat0_inv, mat_inv));
//}
//
//
//// test wrt u2
//BOOST_AUTO_TEST_CASE(L2U2Derivative)
//{
//    DenseRealMatrix eig(5,5), mat(5,5), mat_inv(5,5);
//    RealMatrixX eig0(5,5), mat0(5,5), mat0_inv(5,5);
//    
//    // initialize the analytical derivative for comparison
//    mat0 <<
//    0.5933908290969268, 81.49483049671491, 0.2373563316387707, 776.6885867287626, -754.6885867287625, 0, 0, 0, 0, 0, 0, -1, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
//    mat0_inv <<
//    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.00033233632436025255, -3.0212393123659325e-6, 0.00006646726487205052, 1.5106196561829662e-6, 1.5106196561829662e-6, 0, 0, 0, 0, 0, -0.4176517890908552, -0.00003323363243602526, -0.1677918555499347, -0.00044701814733376236, 0.0004802517797697876;
//    
//    libMesh::Point n(0.5,sqrt(1.-.5*.5-.2*.2),0.2);
//    _fluid_elem->calculate_advection_left_eigenvector_and_inverse_u2_derivative_for_normal(*_primitive_sol,
//                                                                                           n,
//                                                                                           eig,
//                                                                                           mat,
//                                                                                           mat_inv);
//    
//    //BOOST_CHECK(compare(eig0, eig));
//    BOOST_CHECK(compare(mat0, mat));
//    BOOST_CHECK(compare(mat0_inv, mat_inv));
//}
//
//
//// test wrt u3
//BOOST_AUTO_TEST_CASE(L2U3Derivative)
//{
//    DenseRealMatrix eig(5,5), mat(5,5), mat_inv(5,5);
//    RealMatrixX eig0(5,5), mat0(5,5), mat0_inv(5,5);
//    
//    // initialize the analytical derivative for comparison
//    mat0 <<
//    0, -19.389080351973522, -1, 203.74103412361382, -159.74103412361384, 0, 0, 0, 0, 0, 0, -0.2373563316387707, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0;
//    mat0_inv <<
//    0, 0, 3.0212393123659325e-6, 0, 0, 0, 0, 0.00033233632436025255, 0, 0, 0, 0, 0.00003323363243602526, 0, 0, 0.00033233632436025255, -3.0212393123659325e-6, 0.00013293452974410104, 1.5106196561829662e-6, 1.5106196561829662e-6, -0.09268860086407446, -0.00006646726487205052, 0.9806547025589897, -0.00007681306173209923, 0.00014328032660414973;
//    
//    libMesh::Point n(0.5,sqrt(1.-.5*.5-.2*.2),0.2);
//    _fluid_elem->calculate_advection_left_eigenvector_and_inverse_u3_derivative_for_normal(*_primitive_sol,
//                                                                                           n,
//                                                                                           eig,
//                                                                                           mat,
//                                                                                           mat_inv);
//    
//    //BOOST_CHECK(compare(eig0, eig));
//    BOOST_CHECK(compare(mat0, mat));
//    BOOST_CHECK(compare(mat0_inv, mat_inv));
//}
//
//
//// test wrt T
//BOOST_AUTO_TEST_CASE(L2TDerivative)
//{
//    DenseRealMatrix eig(5,5), mat(5,5), mat_inv(5,5);
//    RealMatrixX eig0(5,5), mat0(5,5), mat0_inv(5,5);
//    
//    // initialize the analytical derivative for comparison
//    mat0 <<
//    0, -1003, 0, 94.54494180110561, -94.54494180110561, 0, 0, 0, -0.6884130080439917, 0.6884130080439917, 0, 0, 0, -1.1601342223163071, 1.1601342223163071, 0, 0, 0, -0.2753652032175967, 0.2753652032175967, 0, 0, 0, 0, 0;
//    mat0_inv <<
//    -1.0070797707886442e-6, 9.15527064353313e-9, -2.0141595415772884e-7, -4.577635321766565e-9, -4.577635321766565e-9, -0.00011077877478675086, 1.0070797707886442e-6, -0.000022155754957350172, -8.6696346878699e-8, -9.203834239099449e-7, -0.000011077877478675086, 1.0070797707886442e-7, -2.2155754957350172e-6, 6.521232289632571e-7, -7.528312060421217e-7, -0.000022155754957350172, 2.0141595415772884e-7, -4.4311509914700345e-6, 6.602943832738479e-8, -2.6744539248511365e-7, -0.006397474243934862, 0.0000581588567630442, -0.0012794948487869725, 0.000028168833386663504, -0.0000863276901497077;
//    
//    libMesh::Point n(0.5,sqrt(1.-.5*.5-.2*.2),0.2);
//    _fluid_elem->calculate_advection_left_eigenvector_and_inverse_T_derivative_for_normal(*_primitive_sol,
//                                                                                          n,
//                                                                                          eig,
//                                                                                          mat,
//                                                                                          mat_inv);
//    
//    //BOOST_CHECK(compare(eig0, eig));
//    BOOST_CHECK(compare(mat0, mat));
//    BOOST_CHECK(compare(mat0_inv, mat_inv));
//}
//
//
//BOOST_AUTO_TEST_SUITE_END()
//
//
//
//BOOST_FIXTURE_TEST_SUITE(L3Derivative, BuildFluidElem)
//
//
//// sensitivity wrt rho is all zero. Hence, the first one tested is wrt u1
//BOOST_AUTO_TEST_CASE(L3U1Derivative)
//{
//    DenseRealMatrix eig(5,5), mat(5,5), mat_inv(5,5);
//    RealMatrixX eig0(5,5), mat0(5,5), mat0_inv(5,5);
//    
//    // initialize the analytical derivative for comparison
//    mat0 <<
//    -1, 0, -96.94540175986761, 564.3525853090346, -344.3525853090345, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, -0.5933908290969268, 0, 0, 0, 0, 0, 0, 0;
//    mat0_inv <<
//    3.0212393123659325e-6, 0, 0, 0, 0, 0.0006646726487205051, 0.00003323363243602526, -3.0212393123659325e-6, 1.5106196561829662e-6, 1.5106196561829662e-6, 0.00003323363243602526, 0, 0, 0, 0, 0.00006646726487205052, 0, 0, 0, 0, 0.8057494184114324, -0.09634430043203723, -0.00033233632436025255, -0.00010894857324018491, 0.0004412848976004374;
//    
//    libMesh::Point n(0.5,0.2,sqrt(1.-.5*.5-.2*.2));
//    _fluid_elem->calculate_advection_left_eigenvector_and_inverse_u1_derivative_for_normal(*_primitive_sol,
//                                                                                           n,
//                                                                                           eig,
//                                                                                           mat,
//                                                                                           mat_inv);
//    
//    //BOOST_CHECK(compare(eig0, eig));
//    BOOST_CHECK(compare(mat0, mat));
//    BOOST_CHECK(compare(mat0_inv, mat_inv));
//}
//
//
//// test wrt u2
//BOOST_AUTO_TEST_CASE(L3U2Derivative)
//{
//    DenseRealMatrix eig(5,5), mat(5,5), mat_inv(5,5);
//    RealMatrixX eig0(5,5), mat0(5,5), mat0_inv(5,5);
//    
//    // initialize the analytical derivative for comparison
//    mat0 <<
//    0, -1, -5.778160703947044, 192.74103412361382, -170.74103412361382, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, -0.2373563316387707, 0, 0, 0, 0, 0, 0, 0;
//    mat0_inv <<
//    0, 3.0212393123659325e-6, 0, 0, 0, 0, 0.00033233632436025255, 0, 0, 0, 0.00033233632436025255, 0.00006646726487205052, -3.0212393123659325e-6, 1.5106196561829662e-6, 1.5106196561829662e-6, 0, 0.00006646726487205052, 0, 0, 0, -0.09634430043203723, 0.9795579926886008, -0.00003323363243602526, -0.00009342987795011185, 0.00012666351038613708;
//    
//    libMesh::Point n(0.5,0.2,sqrt(1.-.5*.5-.2*.2));
//    _fluid_elem->calculate_advection_left_eigenvector_and_inverse_u2_derivative_for_normal(*_primitive_sol,
//                                                                                           n,
//                                                                                           eig,
//                                                                                           mat,
//                                                                                           mat_inv);
//    
//    //BOOST_CHECK(compare(eig0, eig));
//    BOOST_CHECK(compare(mat0, mat));
//    BOOST_CHECK(compare(mat0_inv, mat_inv));
//}
//
//
//// test wrt u3
//BOOST_AUTO_TEST_CASE(L3U3Derivative)
//{
//    DenseRealMatrix eig(5,5), mat(5,5), mat_inv(5,5);
//    RealMatrixX eig0(5,5), mat0(5,5), mat0_inv(5,5);
//    
//    // initialize the analytical derivative for comparison
//    mat0 <<
//    0.5933908290969268, 0.2373563316387707, 89.88391084868843, 787.6885867287626, -743.6885867287626, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, -1, 0, 0, 0, 0, 0;
//    mat0_inv <<
//    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.00033233632436025255, 0.00003323363243602526, -3.0212393123659325e-6, 1.5106196561829662e-6, 1.5106196561829662e-6, -0.4139960895228924, -0.16779185554993467, -0.00006646726487205052, -0.00043040133111574974, 0.0004968685959878003;
//    libMesh::Point n(0.5,0.2,sqrt(1.-.5*.5-.2*.2));
//    _fluid_elem->calculate_advection_left_eigenvector_and_inverse_u3_derivative_for_normal(*_primitive_sol,
//                                                                                           n,
//                                                                                           eig,
//                                                                                           mat,
//                                                                                           mat_inv);
//    
//    //BOOST_CHECK(compare(eig0, eig));
//    BOOST_CHECK(compare(mat0, mat));
//    BOOST_CHECK(compare(mat0_inv, mat_inv));
//}
//
//
//// test wrt T
//BOOST_AUTO_TEST_CASE(L3TDerivative)
//{
//    DenseRealMatrix eig(5,5), mat(5,5), mat_inv(5,5);
//    RealMatrixX eig0(5,5), mat0(5,5), mat0_inv(5,5);
//    
//    // initialize the analytical derivative for comparison
//    mat0 <<
//    0, 0, -1003, 104.27740101119142, -104.27740101119142, 0, 0, 0, -0.6884130080439917, 0.6884130080439917, 0, 0, 0, -0.2753652032175967, 0.2753652032175967, 0, 0, 0, -1.1601342223163071, 1.1601342223163071, 0, 0, 0, 0, 0;
//    mat0_inv <<
//    -1.0070797707886442e-6, -1.0070797707886442e-7, 9.15527064353313e-9, -4.577635321766565e-9, -4.577635321766565e-9, -0.00011077877478675086, -0.000011077877478675086, 1.0070797707886442e-6, -8.6696346878699e-8, -9.203834239099449e-7, -0.000011077877478675086, -1.1077877478675086e-6, 1.0070797707886442e-7, 1.16383426866817e-7, -2.1709140394568142e-7, -0.000022155754957350172, -2.2155754957350172e-6, 2.0141595415772884e-7, 6.01769240423825e-7, -8.031851945815539e-7, -0.006397474243934862, -0.0006397474243934863, 0.0000581588567630442, 0.000034061971209724335, -0.00009222082797276853;
//    
//    libMesh::Point n(0.5,0.2,sqrt(1.-.5*.5-.2*.2));
//    _fluid_elem->calculate_advection_left_eigenvector_and_inverse_T_derivative_for_normal(*_primitive_sol,
//                                                                                          n,
//                                                                                          eig,
//                                                                                          mat,
//                                                                                          mat_inv);
//    
//    //BOOST_CHECK(compare(eig0, eig));
//    BOOST_CHECK(compare(mat0, mat));
//    BOOST_CHECK(compare(mat0_inv, mat_inv));
//}
//
//
//BOOST_AUTO_TEST_SUITE_END()
//

