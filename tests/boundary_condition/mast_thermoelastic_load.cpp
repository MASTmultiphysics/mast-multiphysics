#include "catch.hpp"

#include "test_helpers.h"

// MAST includes
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "base/boundary_condition_base.h"

/** 
 * NOTE: In this test, we aren't really checking anything new that wasn't
 * already checked in the function_set_base and boundary_condition_base tests.
 * We are more or least double checking that the strings specific to a 
 * thermoelastic load work.
 */
TEST_CASE("create_thermoelastic_load",
          "[thermoelastic],[load]")
{
    MAST::Parameter temperature("T", 482.222);
    MAST::Parameter ref_temperature("T0", 0.0);
    
    MAST::ConstantFieldFunction temperature_f("temperature", temperature);
    MAST::ConstantFieldFunction ref_temperature_f("ref_temperature", ref_temperature);
    
    MAST::BoundaryConditionBase temperature_load(MAST::TEMPERATURE);   // Initialize a thermoelastic load
    
    REQUIRE( temperature_load.type() == MAST::TEMPERATURE );
    
    temperature_load.add(temperature_f);                               // Add the field function defining the temperature to the thermoelastic load object
    temperature_load.add(ref_temperature_f);                           // Add the field function defining the reference temperature to the thermoelastic load object
    
    REQUIRE( temperature_load.contains("temperature") );
    REQUIRE( temperature_load.contains("ref_temperature") );
    REQUIRE( temperature_load.depends_on(temperature) );
    REQUIRE( temperature_load.depends_on(ref_temperature) );
    
    REQUIRE_NOTHROW( temperature_load.get<MAST::FieldFunction<Real>>("temperature") );
    REQUIRE_NOTHROW( temperature_load.get<MAST::FieldFunction<Real>>("ref_temperature") );
}
