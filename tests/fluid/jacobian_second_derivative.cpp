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
//#define BOOST_TEST_MODULE InviscidFluxSecondDerivative
//#include <boost/test/unit_test.hpp>
//
//// Test includes
//#include "fluid/fluid_elem_initialization.h"
//
//
//
//BOOST_FIXTURE_TEST_SUITE(fxSecondDerivative, BuildFluidElem)
//
//
//// sensitivity wrt rho is all zero. Hence, the first one tested is wrt u1
//BOOST_AUTO_TEST_CASE(AxU1Derivative)
//{
//    DenseRealMatrix deriv(5,5);
//    RealMatrixX deriv_analytical(5,5);
//    
//    // initialize the analytical derivative of the flux Jacobian
//    deriv_analytical
//    << 0,    0,    0,    0,    0,
//    -(62975./358.),    1145./716.,    0,    0,    0,
//    -11.,    0,    1.,    0,    0,
//    -22.,    0,     0,     1.,     0,
//    -(489809925./1432.), 3905./179., -(3157./716.), -(3157./358.), 1003./716.;
//    
//    _fluid_elem->calculate_advection_flux_jacobian_u1_derivative(0,
//                                                                 *_primitive_sol,
//                                                                 deriv);
//    
//    BOOST_CHECK(compare(deriv_analytical, deriv));
//    }
//    
//    
//    // test wrt u2
//    BOOST_AUTO_TEST_CASE(AxU2Derivative)
//    {
//        DenseRealMatrix deriv(5,5);
//        RealMatrixX deriv_analytical(5,5);
//        
//        // initialize the analytical derivative of the flux Jacobian
//        deriv_analytical <<
//        0, 0, 0, 0, 0,
//        4.409217877094972, 0, -0.40083798882681565, 0, 0,
//        -110, 1, 0, 0, 0,
//        0, 0, 0, 0, 0,
//        -724.9860335195531, 11, -44.09217877094972, 0, 0;
//        
//        _fluid_elem->calculate_advection_flux_jacobian_u2_derivative(0,
//                                                                     *_primitive_sol,
//                                                                     deriv);
//        
//        BOOST_CHECK(compare(deriv_analytical, deriv));
//    }
//    
//    
//    // test wrt u3
//    BOOST_AUTO_TEST_CASE(AxU3Derivative)
//    {
//        DenseRealMatrix deriv(5,5);
//        RealMatrixX deriv_analytical(5,5);
//        
//        // initialize the analytical derivative of the flux Jacobian
//        deriv_analytical <<
//        0, 0, 0, 0, 0,
//        8.818435754189943, 0, 0, -0.40083798882681565, 0,
//        0, 0, 0, 0, 0,
//        -110, 1, 0, 0, 0,
//        -1449.9720670391062, 22, 0, -44.09217877094972, 0;
//        
//        _fluid_elem->calculate_advection_flux_jacobian_u3_derivative(0,
//                                                                     *_primitive_sol,
//                                                                     deriv);
//        
//        BOOST_CHECK(compare(deriv_analytical, deriv));
//    }
//    
//    
//    
//    // test wrt u3
//    BOOST_AUTO_TEST_CASE(AxTDerivative)
//    {
//        DenseRealMatrix deriv(5,5);
//        RealMatrixX deriv_analytical(5,5);
//        
//        // initialize the analytical derivative of the flux Jacobian
//        deriv_analytical <<
//        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -110330., 1003., 0, 0, 0;
//        
//        _fluid_elem->calculate_advection_flux_jacobian_T_derivative(0,
//                                                                    *_primitive_sol,
//                                                                    deriv);
//        
//        BOOST_CHECK(compare(deriv_analytical, deriv));
//    }
//    
//    
//    BOOST_AUTO_TEST_SUITE_END()
//    
//    
//    
//    BOOST_FIXTURE_TEST_SUITE(fySecondDerivative, BuildFluidElem)
//    
//    
//    // sensitivity wrt rho is all zero. Hence, the first one tested is wrt u1
//    BOOST_AUTO_TEST_CASE(AyU1Derivative)
//    {
//        DenseRealMatrix deriv(5,5);
//        RealMatrixX deriv_analytical(5,5);
//        
//        // initialize the analytical derivative of the flux Jacobian
//        deriv_analytical <<
//        0, 0, 0, 0, 0, -11, 0, 1, 0, 0, 44.09217877094972, -0.40083798882681565, 0, 0, 0, 0, 0, 0, 0, 0, -724.9860335195531, -4.409217877094972, 110, 0, 0;
//        
//        _fluid_elem->calculate_advection_flux_jacobian_u1_derivative(1,
//                                                                     *_primitive_sol,
//                                                                     deriv);
//        
//        BOOST_CHECK(compare(deriv_analytical, deriv));
//    }
//    
//    
//    // test wrt u2
//    BOOST_AUTO_TEST_CASE(AyU2Derivative)
//    {
//        DenseRealMatrix deriv(5,5);
//        RealMatrixX deriv_analytical(5,5);
//        
//        // initialize the analytical derivative of the flux Jacobian
//        deriv_analytical <<
//        0, 0, 0, 0, 0, -110, 1, 0, 0, 0, -17.59078212290503, 0, 1.5991620111731844, 0, 0, -22, 0, 0, 1, 0, -334868.6752793296, -44.09217877094972, 2.1815642458100557, -8.818435754189943, 1.4008379888268156;
//        
//        _fluid_elem->calculate_advection_flux_jacobian_u2_derivative(1,
//                                                                     *_primitive_sol,
//                                                                     deriv);
//        
//        BOOST_CHECK(compare(deriv_analytical, deriv));
//    }
//    
//    
//    // test wrt u3
//    BOOST_AUTO_TEST_CASE(AyU3Derivative)
//    {
//        DenseRealMatrix deriv(5,5);
//        RealMatrixX deriv_analytical(5,5);
//        
//        // initialize the analytical derivative of the flux Jacobian
//        deriv_analytical <<
//        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8.818435754189943, 0, 0, -0.40083798882681565, 0, -11, 0, 1, 0, 0, -144.99720670391062, 0, 22, -4.409217877094972, 0;
//        
//        _fluid_elem->calculate_advection_flux_jacobian_u3_derivative(1,
//                                                                     *_primitive_sol,
//                                                                     deriv);
//        
//        BOOST_CHECK(compare(deriv_analytical, deriv));
//    }
//    
//    
//    // test wrt T
//    BOOST_AUTO_TEST_CASE(AyTDerivative)
//    {
//        DenseRealMatrix deriv(5,5);
//        RealMatrixX deriv_analytical(5,5);
//        
//        // initialize the analytical derivative of the flux Jacobian
//        deriv_analytical <<
//        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -11033, 0, 1003, 0, 0;
//        
//        _fluid_elem->calculate_advection_flux_jacobian_T_derivative(1,
//                                                                    *_primitive_sol,
//                                                                    deriv);
//        
//        BOOST_CHECK(compare(deriv_analytical, deriv));
//    }
//    
//    
//    BOOST_AUTO_TEST_SUITE_END()
//
//    
//    
//    BOOST_FIXTURE_TEST_SUITE(fzSecondDerivative, BuildFluidElem)
//    
//    
//    // sensitivity wrt rho is all zero. Hence, the first one tested is wrt u1
//    BOOST_AUTO_TEST_CASE(AzU1Derivative)
//    {
//        DenseRealMatrix deriv(5,5);
//        RealMatrixX deriv_analytical(5,5);
//        
//        // initialize the analytical derivative of the flux Jacobian
//        deriv_analytical <<
//        0, 0, 0, 0, 0, -22, 0, 0, 1, 0, 0, 0, 0, 0, 0, 44.09217877094972, -0.40083798882681565, 0, 0, 0, -1449.9720670391062, -8.818435754189943, 0, 110, 0;
//        
//        _fluid_elem->calculate_advection_flux_jacobian_u1_derivative(2,
//                                                                     *_primitive_sol,
//                                                                     deriv);
//        
//        BOOST_CHECK(compare(deriv_analytical, deriv));
//    }
//    
//    
//    // test wrt u2
//    BOOST_AUTO_TEST_CASE(AzU2Derivative)
//    {
//        DenseRealMatrix deriv(5,5);
//        RealMatrixX deriv_analytical(5,5);
//        
//        // initialize the analytical derivative of the flux Jacobian
//        deriv_analytical <<
//        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -22, 0, 0, 1, 0, 4.409217877094972, 0, -0.40083798882681565, 0, 0, -144.99720670391062, 0, -8.818435754189943, 11, 0;
//        
//        _fluid_elem->calculate_advection_flux_jacobian_u2_derivative(2,
//                                                                     *_primitive_sol,
//                                                                     deriv);
//        
//        BOOST_CHECK(compare(deriv_analytical, deriv));
//    }
//    
//    
//    // test wrt u3
//    BOOST_AUTO_TEST_CASE(AzU3Derivative)
//    {
//        DenseRealMatrix deriv(5,5);
//        RealMatrixX deriv_analytical(5,5);
//        
//        // initialize the analytical derivative of the flux Jacobian
//        deriv_analytical <<
//        0, 0, 0, 0, 0, -110, 1, 0, 0, 0, -11, 0, 1, 0, 0, -35.18156424581006, 0, 0, 1.5991620111731844, 0, -335086.17108938546, -44.09217877094972, -4.409217877094972, 4.363128491620111, 1.4008379888268156;
//        
//        _fluid_elem->calculate_advection_flux_jacobian_u3_derivative(2,
//                                                                     *_primitive_sol,
//                                                                     deriv);
//        
//        BOOST_CHECK(compare(deriv_analytical, deriv));
//    }
//    
//    
//    // test wrt T
//    BOOST_AUTO_TEST_CASE(AzTDerivative)
//    {
//        DenseRealMatrix deriv(5,5);
//        RealMatrixX deriv_analytical(5,5);
//        
//        // initialize the analytical derivative of the flux Jacobian
//        deriv_analytical <<
//        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -22066, 0, 0, 1003, 0;
//        
//        _fluid_elem->calculate_advection_flux_jacobian_T_derivative(2,
//                                                                    *_primitive_sol,
//                                                                    deriv);
//        
//        BOOST_CHECK(compare(deriv_analytical, deriv));
//    }
//    
//    
//    BOOST_AUTO_TEST_SUITE_END()
//
//
