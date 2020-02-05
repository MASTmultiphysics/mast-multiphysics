// Catch2 includes
#include "catch.hpp"

// libMesh includes
#include "libmesh/point.h"

// MAST includes
#include "base/parameter.h"
#include "base/constant_field_function.h"

TEST_CASE("constant_field_functions", "[field_function],[constant],[base]")
{
    /**
     * This code, outside of any 'SECTION' is run automatically each time before
     * a section is run.
     */
    
    // Create a parameter
    const Real initial_value = 4.984;
    MAST::Parameter parameter1("p1", initial_value);
    
    const std::string initial_name = "cf1";
    MAST::ConstantFieldFunction cfield1(initial_name, parameter1);
    
    /**
     * Only one section is run at a time. Each time a section completes, the 
     * test starts over at the beginning of the 'TEST_CASE' before running
     * the next section.
     */
    
    SECTION("constant field function can return a constant reference to its name")
    {
        const std::string& name = cfield1.name();
        CHECK( name == initial_name );
    }
    
    SECTION("constant field function can return a copy of its name")
    {
        std::string name = cfield1.name();
        CHECK( name == initial_name );
        name += "_added_string";
        CHECK( name != initial_name );
        CHECK ( cfield1.name() == initial_name );
    }
    
    SECTION("constant field function can return the value")
    {
        Real cfield_value;
        cfield1(cfield_value);
        CHECK( cfield_value == initial_value );
    }
    
    SECTION("constant field function's derivative w.r.t. its own parameter is 1")
    {
        Real dvalue_dparam;
        cfield1.derivative(parameter1, dvalue_dparam);
        REQUIRE ( dvalue_dparam == 1.0 );
    }
    
    SECTION("constant field function's derivative w.r.t. another parameter is 0")
    {
        Real dvalue_dparam;
        MAST::Parameter parameter2("p2", 3.623);
        cfield1.derivative(parameter2, dvalue_dparam);
        CHECK( dvalue_dparam == 0.0 );
    }
    
    SECTION("constant field function can return the value given a point and time")
    {
        Real cfield_value;
        const libMesh::Point point(2.3, 3.1, 5.2);
        const Real time = 2.34;
        cfield1(point, time, cfield_value);
        CHECK( cfield_value == initial_value );
    }
    
    SECTION("constant field function is constant over time and space")
    {
        Real cfield_value;
        const libMesh::Point point(2.3, 3.1, 5.2);
        const Real time = 2.34;
        cfield1(point, time, cfield_value);
        CHECK( cfield_value == initial_value );
        
        Real cfield_value2;
        const libMesh::Point point2(4.3, -16.4, 55.2);
        const Real time2 = 88.4;
        cfield1(point2, time2, cfield_value2);
        CHECK( cfield_value2 == initial_value );
        
        Real cfield_value3;
        const libMesh::Point point3(-18.0, 10.2, -7.5);
        const Real time3 = 36.7;
        cfield1(point3, time3, cfield_value3);
        CHECK( cfield_value3 == initial_value );
    }
    
    SECTION("constant field function can return derivative to its own parameter at a given point and time")
    {
        Real dvalue_dparam;
        const libMesh::Point point(3.5, -2.4, -11.2);
        const Real time = 6.87;
        cfield1.derivative(parameter1, point, time, dvalue_dparam);
        CHECK( dvalue_dparam == 1.0 );
    }
    
    SECTION("constant field function can return derivative to other parameter at a given point and time")
    {
        MAST::Parameter parameter2("p2", 5.478);
        Real dvalue_dparam;
        const libMesh::Point point(3.5, -2.4, -11.2);
        const Real time = 6.87;
        cfield1.derivative(parameter2, point, time, dvalue_dparam);
        CHECK( dvalue_dparam == 0.0 );
    }
    
    
    SECTION("constant field function derivative w.r.t. itself is constant 1 over time and space")
    {
        Real dvalue_dparam;
        const libMesh::Point point(-5.0, 8.4, 7.3);
        const Real time = 3.5;
        cfield1.derivative(parameter1, point, time, dvalue_dparam);
        CHECK( dvalue_dparam == 1.0 );
        
        Real dvalue_dparam2;
        const libMesh::Point point2(8.8, -10.5, 200.3);
        const Real time2 = 107.5;
        cfield1.derivative(parameter1, point2, time2, dvalue_dparam2);
        CHECK( dvalue_dparam2 == 1.0 );
        
        Real dvalue_dparam3;
        const libMesh::Point point3(0.0, 57.8,-150.7);
        const Real time3 = 0.0;
        cfield1.derivative(parameter1, point3, time3, dvalue_dparam3);
        CHECK( dvalue_dparam3 == 1.0 );
    }
    
    SECTION("constant field function derivative w.r.t. other parameters is constant 0 over time and space")
    {
        MAST::Parameter parameter2("p2", 6.578);
        
        Real dvalue_dparam;
        const libMesh::Point point(-5.0, 8.4, 7.3);
        const Real time = 3.5;
        cfield1.derivative(parameter2, point, time, dvalue_dparam);
        CHECK( dvalue_dparam == 0.0 );
        
        Real dvalue_dparam2;
        const libMesh::Point point2(8.8, -10.5, 200.3);
        const Real time2 = 107.5;
        cfield1.derivative(parameter2, point2, time2, dvalue_dparam2);
        CHECK( dvalue_dparam2 == 0.0 );
        
        Real dvalue_dparam3;
        const libMesh::Point point3(0.0, 57.8,-150.7);
        const Real time3 = 0.0;
        cfield1.derivative(parameter2, point3, time3, dvalue_dparam3);
        CHECK( dvalue_dparam3 == 0.0 );
    }
    
    SECTION("constant field function is NOT a shape parameter by default")
    {
        CHECK_FALSE( cfield1.is_shape_parameter() );
    }
    
    SECTION("constant field function can be set as a shape parameter")
    {
        cfield1.set_as_shape_parameter(true);
        CHECK( cfield1.is_shape_parameter() );
    }
    
    SECTION("constant field function is NOT a topology parameter by default")
    {
        CHECK_FALSE( cfield1.is_topology_parameter() );
    }
    
    SECTION("constant field function can be set as a topology parameter")
    {
        cfield1.set_as_topology_parameter(true);
        CHECK( cfield1.is_topology_parameter() );
    }
}
