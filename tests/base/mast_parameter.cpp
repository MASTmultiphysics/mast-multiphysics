// Catch2 includes
#include "catch.hpp"

// MAST includes
#include "base/parameter.h"


/**
 * MAST::Parameter objects are independent of libMesh and thus we do not need
 * to initialize any libMesh or MPI communicator objects here.
 */
TEST_CASE("parameters", "[parameter],[base]")
{
    /**
     * This code, outside of any 'SECTION' is run automatically each time before
     * a section is run.
     */
    
    // Create a parameter
    const Real initial_value = 4.984;
    const std::string initial_name = "p1";
    MAST::Parameter parameter1(initial_name, initial_value);
    
    
    /**
     * Only one section is run at a time. Each time a section completes, the 
     * test starts over at the beginning of the 'TEST_CASE' before running
     * the next section.
     */
    SECTION("parameter can return a constant reference to its name")
    {
        const std::string& name = parameter1.name();
        CHECK( name == initial_name );
    }
    
    SECTION("parameter can return a copy of its name")
    {
        std::string name = parameter1.name();
        CHECK( name == initial_name );
        name += "_added_string";
        CHECK( name != initial_name );
        CHECK ( parameter1.name() == initial_name );
    }
    
    SECTION("parameter can return a constant reference to its value")
    {
        const Real const_param_value = parameter1();
        CHECK( const_param_value == initial_value );
    }
    
    SECTION("parameter can return a writable reference to its value")
    {
        Real& parameter_value = parameter1();
        CHECK( parameter_value == initial_value );
    }
    
    SECTION("parameter value can be changed through writable reference to its value")
    {
        Real& parameter_value = parameter1();
        parameter_value = 2.4578;
        CHECK( parameter1() == 2.4578 );
    }
    
    SECTION("parameter value can be set through assignment oeprator '='")
    {
        parameter1 = 5.678;
        CHECK( parameter1() == 5.678 );
    }
    
    SECTION("parameter can return a pointer to its value")
    {
        Real* val_ptr = parameter1.ptr();
        CHECK( *val_ptr == initial_value );
    }
    
    SECTION("parameter depends on itself")
    {
        const bool depends_on_itself = parameter1.depends_on(parameter1);
        CHECK( depends_on_itself );
    }
    
    SECTION("parameter does not depend on other field functions")
    {
        MAST::Parameter parameter2("p2", 5.242);
        const bool depends_on_other = parameter1.depends_on(parameter2);
        CHECK_FALSE( depends_on_other );
    }
    
    SECTION("parameter is NOT a shape parameter by default")
    {
        CHECK_FALSE( parameter1.is_shape_parameter() );
    }
    
    SECTION("parameter can be set as a shape parameter")
    {
        parameter1.set_as_shape_parameter(true);
        CHECK( parameter1.is_shape_parameter() );
    }
    
    SECTION("parameter is NOT a topology parameter by default")
    {
        CHECK_FALSE( parameter1.is_topology_parameter() );
    }
    
    SECTION("parameter can be set as a topology parameter")
    {
        parameter1.set_as_topology_parameter(true);
        CHECK( parameter1.is_topology_parameter() );
    }
    
    //TODO: Test that Parameter constructor accepts a field function as well
}
