project (MAST)
enable_language (CXX Fortran)

set (MAST_VERSION_MAJOR 0)
set (MAST_VERSION_MINOR 9)

cmake_minimum_required(VERSION 2.8)

# flags for optional compilation of packages
option (ENABLE_GCMMA "Build with GCMMA interface" OFF)
option (ENABLE_DOT   "Build with DOT interface"   OFF)
option (ENABLE_NPSOL  "Build with NPSOL interface" OFF)
option (ENABLE_MATPLOTLIB_CPP  "Build with Matplotlib-cpp interface" OFF)
option (ENABLE_CYTHON  "Build with CYTHON interface" OFF)

#set default values
set (MAST_ENABLE_GCMMA 0)
set (MAST_ENABLE_DOT   0)
set (MAST_ENABLE_NPSOL 0)
set (MAST_ENABLE_MATPLOTLIB   0)
set (MAST_ENABLE_CYTHON 0)
set (MAST_ENABLE_PYTHON 0)

# look for gcmma library

if (ENABLE_GCMMA)
   set (MAST_ENABLE_GCMMA 1)
   set (gcmma_lib_file    "gcmma-lib-file" CACHE STRING "GCMMA optimization library")
   set (gcmma_lib_dir     "gcmma-lib-dir" CACHE STRING "Directory containing GCMMA optimization library")
   find_library (gcmma_lib   ${gcmma_lib_file}    ${gcmma_lib_dir})
endif()


# look for dot library
if (ENABLE_DOT)
    set (MAST_ENABLE_DOT 1)
    set (dot_lib_file      "dot-lib-file" CACHE STRING "DOT optimization library")
    set (dot_lib_dir       "dot-lib-dir" CACHE STRING "Directory containing DOT optimization library")
    find_library (dot_lib     ${dot_lib_file}      ${dot_lib_dir})
endif()


# look for npsol library
if (ENABLE_NPSOL)
   set (MAST_ENABLE_NPSOL 1)
   set (npsol_lib_file    "npsol-lib-file" CACHE STRING "NPSOL optimization library")
   set (npsol_lib_dir     "npsol-lib-dir" CACHE STRING "Directory containing NPSOL optimization library")
   find_library (npsol_lib   ${npsol_lib_file}    ${npsol_lib_dir})
endif()


#look for cython and python
if (ENABLE_CYTHON)
   set (MAST_ENABLE_CYTHON 1)
   set (MAST_ENABLE_PYTHON 1)
   set (CYTHON            "cython-compiler" CACHE STRING "Cython executable")
endif()


#look for matplotlib files 
if (ENABLE_MATPLOTLIB_CPP)
   set (MAST_ENABLE_MATPLOTLIB 1)
   set (MAST_ENABLE_PYTHON 1)
   set (matplotlib_cpp_include_dir "matplot-cpp-include-dir" CACHE STRING "Matplotlib-cpp include directory") 
endif()

if (MAST_ENABLE_PYTHON) 
   set (python_include_dir  "python-include-dir" CACHE STRING "Directory containing python include files")
   set (python_lib_file     "python-lib-file" CACHE STRING "Python library to link against")
   set (python_lib_dir      "python-lib-dir" CACHE STRING "Directory containing python library")
   find_library (python_lib   ${python_lib_file}    ${python_lib_dir} NO_DEFAULT_PATH)
endif()


set (libmesh_dir       "libmesh-dir" CACHE STRING "Directory containing libMesh header and library")
set (eigen_include_dir "eigen-include-dir" CACHE STRING "Directory containing eigen header files")
set (boost_include_dir "boost-include-dir" CACHE STRING "Directory containing boost header files")
set (boost_test_lib    "boost-unit-test-lib" CACHE STRING "BOOST unit test library")
set (boost_lib_dir     "boost-lib-dir" CACHE STRING "Directory containing boost unit test library files")
set (mpi_include_dir   "mpi-include-dir" CACHE STRING "Directory containing MPI header files")
set (mpi_lib_dir       "mpi-lib-dir" CACHE STRING "Directory containing MPI library files")
set (mpicxx_lib_file   "mpi_cxx" CACHE STRING "Root filename of MPI C++ library. eg. if file is libmpicxx.so, this string is mpicxx")
set (numpy_dir         "numpy-dir" CACHE STRING "Directory containing NUMPY header files and library")
set (fortran_lib_file  "fortran-lib-file" CACHE STRING "fortran library")
set (fortran_lib_dir   "fortran-lib-dir" CACHE STRING "Directory containing fortran library")
set (lapack_lib_file   "lapack" CACHE STRING "LAPACK library")
set (lapack_lib_dir    "lapack-lib-dir" CACHE STRING "Directory containing LAPACK library")
separate_arguments (lapack_lib_file)
set (blas_lib_file     "blas" CACHE STRING "BLAS library")
set (blas_lib_dir      "blas-lib-dir" CACHE STRING "Directory containing BLAS library")
separate_arguments (blas_lib_file)
set (petsc_dir         "petsc-dir" CACHE STRING "Directory containing PETSc installation")
set (slepc_dir         "slepc-dir" CACHE STRING "Directory containing PETSc installation")
set (CMAKE_Fortran_COMPILER "fortran-compiler" CACHE STRING "Fortran compiler")


include_directories (${PROJECT_SOURCE_DIR}/../src
                     ${PROJECT_SOURCE_DIR}/../examples
         		     ${PROJECT_SOURCE_DIR}/../tests
                     ${PROJECT_SOURCE_DIR}/../
                             ${eigen_include_dir}
		             ${boost_include_dir})

include_directories (SYSTEM ${mpi_include_dir})
include_directories (SYSTEM ${numpy_dir}/include)


find_library (mesh_dbg mesh_dbg ${libmesh_dir}/lib)
find_library (mesh_opt mesh_opt ${libmesh_dir}/lib)

find_library (petsc_lib petsc ${petsc_dir}/lib)
find_library (slepc_lib slepc ${slepc_dir}/lib)

find_library (mpicxx_lib    ${mpicxx_lib_file}   ${mpi_lib_dir})
find_library (mpi_lib       mpi                  ${mpi_lib_dir})

find_library (boost_unit_test_lib   ${boost_test_lib}    ${boost_lib_dir})

find_library (fortran_lib   ${fortran_lib_file}    ${fortran_lib_dir})

foreach (f ${lapack_lib_file})
   find_library (lapack_lib    ${f}   ${lapack_lib_dir})
endforeach (f)

foreach (f ${blas_lib_file})
   find_library (blas_lib      ${f}     ${blas_lib_dir})
endforeach (f)


#
# write the configuration file
#
configure_file ( ${PROJECT_SOURCE_DIR}/../src/base/mast_config.h.in
                 ${PROJECT_SOURCE_DIR}/../src/base/mast_config.h )


####################################################################
#  add files for MAST library and add it to cmake's build list
####################################################################
file (GLOB_RECURSE mast_source_files
      ${PROJECT_SOURCE_DIR}/../src/*.cpp
      ${PROJECT_SOURCE_DIR}/../src/*.h)
add_library (mast ${mast_source_files})

target_link_libraries (mast  
		               debug ${mesh_dbg}
		               optimized ${mesh_opt}
		               ${petsc_lib}
		               ${slepc_lib}
		               ${lapack_lib}
                       ${blas_lib}
		               ${gcmma_lib}
		               ${dot_lib}
		               ${npsol_lib}
		               ${fortran_lib}
		               ${mpicxx_lib}
		               ${mpi_lib})

set_property (TARGET mast APPEND 
              PROPERTY INCLUDE_DIRECTORIES
	          ${libmesh_dir}/include)
set_property (TARGET mast APPEND 
	          PROPERTY INCLUDE_DIRECTORIES
	          ${petsc_dir}/include)
set_property (TARGET mast APPEND 
	          PROPERTY INCLUDE_DIRECTORIES
	          ${slepc_dir}/include)
#set_property (TARGET mast APPEND
#              PROPERTY INCLUDE_DIRECTORIES
#              ${numpy_dir}/include)

if (MAST_ENABLE_MATPLOTLIB)
set_property (TARGET mast APPEND
	          PROPERTY INCLUDE_DIRECTORIES
	          ${matplotlib_cpp_include_dir})
endif()

if (MAST_ENABLE_PYTHON)
set_property (TARGET mast APPEND
                PROPERTY INCLUDE_DIRECTORIES
                ${python_include_dir})
target_link_libraries (mast
                        ${python_lib})
endif ()

####################################################################
#  tell cmake to link example driver with mast library
####################################################################

#add the cython compilation
if (ENABLE_CYTHON)
   add_custom_command( OUTPUT 
                     ${PROJECT_SOURCE_DIR}/../examples/structural/nastran_model_analysis/pynastran.h 
  		     ${PROJECT_SOURCE_DIR}/../examples/structural/nastran_model_analysis/pynastran.cpp
                     COMMAND ${CYTHON} --capi-reexport-cincludes --cplus ${PROJECT_SOURCE_DIR}/../examples/structural/nastran_model_analysis/pynastran.pyx  -o ${PROJECT_SOURCE_DIR}/../examples/structural/nastran_model_analysis/pynastran.cpp
                     DEPENDS 
		     ${PROJECT_SOURCE_DIR}/../examples/structural/nastran_model_analysis/pynastran.pyx 
		     ${PROJECT_SOURCE_DIR}/../examples/structural/nastran_model_analysis/mast_interface.pxd
		     IMPLICIT_DEPENDS
		     CXX ${PROJECT_SOURCE_DIR}/../examples/structural/nastran_model_analysis/mast_interface.h
		     CXX ${PROJECT_SOURCE_DIR}/../examples/structural/nastran_model_analysis/mast_interface.cpp 
              	     )    

   set_source_files_properties(${PROJECT_SOURCE_DIR}/../examples/structural/nastran_model_analysis/pynastran.h
			       ${PROJECT_SOURCE_DIR}/../examples/structural/nastran_model_analysis/pynastran.cpp
                               PROPERTIES GENERATED TRUE ) 

   set (pynastran_cpp_file ${PROJECT_SOURCE_DIR}/../examples/structural/nastran_model_analysis/pynastran.cpp)
endif()


file (GLOB_RECURSE example_driver_source_files
      ${PROJECT_SOURCE_DIR}/../examples/*.cpp
      ${PROJECT_SOURCE_DIR}/../examples/*.h)
add_executable (example_driver ${example_driver_source_files} ${pynastran_cpp_file})

target_link_libraries (example_driver
                        mast)
set_property (TARGET example_driver APPEND
              PROPERTY INCLUDE_DIRECTORIES
	          ${libmesh_dir}/include)
set_property (TARGET example_driver APPEND 
              PROPERTY INCLUDE_DIRECTORIES
              ${petsc_dir}/include)
set_property (TARGET example_driver APPEND 
              PROPERTY INCLUDE_DIRECTORIES
              ${slepc_dir}/include)
set_property (TARGET example_driver APPEND
	          PROPERTY INCLUDE_DIRECTORIES
	          ${matplotlib_cpp_include_dir})

if (MAST_ENABLE_PYTHON)
set_property (TARGET example_driver APPEND
                PROPERTY INCLUDE_DIRECTORIES
                ${python_include_dir})
target_link_libraries (example_driver
                    mast
                    ${python_lib})
endif ()

####################################################################
#  tell cmake to link the tests
####################################################################
file (GLOB_RECURSE mast_test_source_files
      ${PROJECT_SOURCE_DIR}/../examples/*.cpp
      ${PROJECT_SOURCE_DIR}/../examples/*.h
      ${PROJECT_SOURCE_DIR}/../tests/*.cpp
      ${PROJECT_SOURCE_DIR}/../tests/*.h)
list (REMOVE_ITEM mast_test_source_files ${PROJECT_SOURCE_DIR}/../examples/base/examples_driver.cpp)
add_executable (mast_tests ${mast_test_source_files})

target_link_libraries (mast_tests
                       mast
			           ${boost_unit_test_lib})
set_property (TARGET mast_tests APPEND 
	          PROPERTY INCLUDE_DIRECTORIES
	          ${libmesh_dir}/include
	          ${})
set_property (TARGET mast_tests APPEND 
              PROPERTY INCLUDE_DIRECTORIES
              ${petsc_dir}/include)
set_property (TARGET mast_tests APPEND 
              PROPERTY INCLUDE_DIRECTORIES
              ${slepc_dir}/include)


