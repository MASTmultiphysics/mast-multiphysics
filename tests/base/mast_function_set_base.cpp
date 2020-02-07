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
