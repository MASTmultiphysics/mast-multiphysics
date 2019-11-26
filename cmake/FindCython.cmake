# Find the Cython compiler.
#
# This code sets the following variables:
#
#  Cython_FOUND - Cython is available on the system.
#  Cython_EXECUTABLE - Path to the Cython executable.
#
# Note this file is a modified version of that provided by Kitware online.
# --- see original Kitware copyright statement at bottom of file.

# Use the Cython executable that lives next to the Python executable if it
# exists. If not, we will use whatever we can find on the system paths.
if(Python3_FOUND)
    get_filename_component(_PYTHON_PATH ${Python3_EXECUTABLE} PATH)
    find_program(Cython_EXECUTABLE
            NAMES
            cython
            cython.bat
            cython3
            HINTS ${_PYTHON_PATH})
else()
    find_program(Cython_EXECUTABLE
            NAMES
            cython
            cython.bat
            cython3)
endif()

# Set standard CMake variables and output status.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Cython
        REQUIRED_VARS
        Cython_EXECUTABLE)

# Advance cache variables.
mark_as_advanced(Cython_EXECUTABLE)





# Original Kitware copywrite statement. Modifications to code included
# formatting changes to utilize MAST-Interface version of FindPython.cmake.
#=============================================================================
# Copyright 2011 Kitware, Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#=============================================================================