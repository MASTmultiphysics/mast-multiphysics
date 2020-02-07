#include "catch.hpp"

#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "base/function_set_base.h"

TEST_CASE("function_set_base",
          "[base][function_set_base]")
{
    MAST::Parameter param("dummy_param", 482.222);
    MAST::Parameter param2("unused_param", 5.5);
    
    MAST::ConstantFieldFunction dummy_f("dummy", param);
    
    MAST::FunctionSetBase f_set;
    
    // Ensure field functions can be added
    REQUIRE_NOTHROW( f_set.add(dummy_f) );
    
    SECTION("function_set_contains")
    {
        // Ensure the correct boolean value is return when checking if a function
        // set depends on a particular field function
        REQUIRE( f_set.contains("dummy") );
        REQUIRE_FALSE( f_set.contains("dummy2") );
    }
    
    SECTION("function_set_depends_on")
    {
        // Ensure the correct boolean value is returned when checking if a function
        // set depends on a particular parameter
        REQUIRE( f_set.depends_on(param) );
        REQUIRE_FALSE( f_set.depends_on(param2) );
    }
    
    SECTION("function_set_get_existing_field_function")
    {
        // Ensure no error is thrown when retrieving a field function from the
        // function set
        REQUIRE_NOTHROW( f_set.get<MAST::FieldFunction<Real>>("dummy") );
    }
    
    SECTION("function_set_get_nonexistant_field_function_error")
    {
        // Ensure an error is thrown when the requested field function doesn't
        // exist in the function set
        REQUIRE_THROWS( f_set.get<MAST::FieldFunction<Real>>("does_not_exist") );
    }
}
