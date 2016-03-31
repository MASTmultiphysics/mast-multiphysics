project (MAST)

cmake_minimum_required(VERSION 2.8)

enable_language (CXX Fortran)

set (MAST_VERSION_MAJOR 0)
set (MAST_VERSION_MINOR 9)

set (libmesh_dir       "libmesh-dir" CACHE STRING "Directory containing libMesh header and library")
#set (vtk_include_dir   "vtk-include-dir" CACHE STRING "Directory containing VTK header files")
#set (vtk_lib_dir       "vtk-lib-dir" CACHE STRING "Directory containing  VTK library files" )
#set (vtk_lib_file      "vtkIO vtkCommon vtkFiltering" CACHE STRING "VTK library files" )
#separate_arguments (vtk_lib_file)
#set (tbb_include_dir   "tbb-include-dir" CACHE STRING "Directory containing TBB header files")
#set (tbb_lib_dir       "tbb-lib-dir" CACHE STRING "Directory containing TBB library files")
#set (tbb_lib_file      "tbb tbbmalloc" CACHE STRING "TBB library files")
#separate_arguments (tbb_lib_file)
set (boost_include_dir "boost-include-dir" CACHE STRING "Directory containing boost header files")
set (boost_test_lib    "boost-unit-test-lib" CACHE STRING "BOOST unit test library")
set (boost_lib_dir     "boost-lib-dir" CACHE STRING "Directory containing boost unit test library files")
set (mpi_include_dir   "mpi-include-dir" CACHE STRING "Directory containing MPI header files")
set (mpi_lib_dir       "mpi-lib-dir" CACHE STRING "Directory containing MPI library files")
set (gcmma_lib_file    "gcmma-lib-file" CACHE STRING "GCMMA optimization library")
set (gcmma_lib_dir     "gcmma-lib-dir" CACHE STRING "Directory containing GCMMA optimization library")
set (dot_lib_file      "dot-lib-file" CACHE STRING "DOT optimization library")
set (dot_lib_dir       "dot-lib-dir" CACHE STRING "Directory containing DOT optimization library")
set (npsol_lib_file    "npsol-lib-file" CACHE STRING "NPSOL optimization library")
set (npsol_lib_dir     "npsol-lib-dir" CACHE STRING "Directory containing NPSOL optimization library")
set (fortran_lib_file  "fortran-lib-file" CACHE STRING "fortran library")
set (fortran_lib_dir   "fortran-lib-dir" CACHE STRING "Directory containing fortran library")
set (lapack_lib_file   "lapack" CACHE STRING "LAPACK library")
separate_arguments (lapack_lib_file)
set (blas_lib_file     "blas" CACHE STRING "BLAS library")
separate_arguments (blas_lib_file)
set (lapack_lib_dir    "lapack-lib-dir" CACHE STRING "Directory containing LAPACK library")
set (blas_lib_dir      "blas-lib-dir" CACHE STRING "Directory containing BLAS library")
set (petsc_dir         "petsc-dir" CACHE STRING "Directory containing PETSc installation")
set (slepc_dir         "slepc-dir" CACHE STRING "Directory containing PETSc installation")
set (CMAKE_Fortran_COMPILER "fortran-compiler" CACHE STRING "Fortran compiler")

#set (EXTRA_LIBS ${EXTRA_LIBS} MathFunctions)


include_directories (../src 
                     ../examples
		     ../tests
		     ../
		    #${vtk_include_dir}
		    #${tbb_include_dir}
		    ${boost_include_dir})

include_directories (SYSTEM ${mpi_include_dir})


find_library (mesh_dbg mesh_dbg ${libmesh_dir}/lib)
find_library (mesh_opt mesh_opt ${libmesh_dir}/lib)

find_library (petsc_lib petsc ${petsc_dir}/lib)
find_library (slepc_lib slepc ${slepc_dir}/lib)

find_library (mpicxx_lib    mpi_cxx              ${mpi_lib_dir})
find_library (mpi_lib       mpi                  ${mpi_lib_dir})

find_library (boost_unit_test_lib   ${boost_test_lib}    ${boost_lib_dir})

find_library (gcmma_lib   ${gcmma_lib_file}    ${gcmma_lib_dir})
find_library (dot_lib     ${dot_lib_file}      ${dot_lib_dir})
find_library (npsol_lib   ${npsol_lib_file}    ${npsol_lib_dir})

find_library (fortran_lib   ${fortran_lib_file}    ${fortran_lib_dir})

#foreach (f ${tbb_lib_file})
#   find_library (tbb_lib    ${f}   ${tbb_lib_dir})
#endforeach (f)

#foreach (f ${vtk_lib_file})
#   find_library (vtk_lib    ${f}   ${vtk_lib_dir})
#endforeach (f)

foreach (f ${lapack_lib_file})
   find_library (lapack_lib    ${f}   ${lapack_lib_dir})
endforeach (f)

foreach (f ${blas_lib_file})
find_library (blas_lib      ${f}     ${blas_lib_dir})
endforeach (f)


#
#  add files for MAST library and add it to cmake's build list
#
file (GLOB_RECURSE mast_source_files ../src/*.cpp ../src/*.h)
add_library (mast ${mast_source_files})

target_link_libraries (mast  
		      debug ${mesh_dbg}
		      optimized ${mesh_opt}
		      ${petsc_lib}
		      ${slepc_lib}
		      ${lapack_lib}
		      ${blas_lib}
		      #${vtk_lib}
		      #${tbb_lib}
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

#
#  tell cmake to link example driver with mast library
#
file (GLOB_RECURSE example_driver_source_files
     	../examples/thermal/*.cpp
	../examples/thermal/*.h
	../examples/fluid/*.cpp 
	../examples/fluid/*.h
	../examples/fsi/*.cpp
	../examples/fsi/*.h
	../examples/structural/*.cpp
	../examples/structural/*.h
	../examples/base/*.cpp) 
add_executable (example_driver ${example_driver_source_files})

target_link_libraries (example_driver
			mast)
set_property (TARGET example_driver APPEND 
	     PROPERTY INCLUDE_DIRECTORIES 
	     ${libmesh_dir}/include) 

#
#  tell cmake to link the tests
#
file (GLOB_RECURSE mast_test_source_files 
     	../examples/thermal/*.cpp
	../examples/thermal/*.h
	../examples/fluid/*.cpp 
	../examples/fluid/*.h
	../examples/base/multilinear_interpolation.cpp
	../examples/structural/bar_extension/*.cpp
	../examples/structural/beam_bending/*.cpp
	../examples/structural/beam_optimization/beam_optimization_base.cpp
	../examples/structural/beam_bending_thermal_stress_with_offset/*.cpp
	../examples/structural/beam_bending_with_offset/*.cpp
	../examples/structural/beam_modal_analysis/*.cpp
	../examples/structural/membrane_extension_biaxial_stress/*.cpp
	../examples/structural/membrane_extension_uniaxial_stress/*.cpp
	../examples/structural/nonlinear_circular_cantilever/*.cpp
	../examples/structural/plate_bending/*.cpp
	../examples/structural/plate_optimization/plate_optimization_base.cpp
	../examples/structural/plate_bending_section_offset/*.cpp
	../examples/structural/plate_bending_thermal_stress/*.cpp
	../examples/structural/plate_modal_analysis/*.cpp
	../examples/structural/plate_thermally_stressed_modal_analysis/*.cpp
	../examples/structural/beam_buckling_prestress/*.cpp
	../examples/structural/beam_piston_theory_flutter/*.cpp
	../examples/structural/plate_piston_theory_flutter/*.cpp
	../examples/structural/plate_thermally_stressed_piston_theory_flutter/*.cpp
	../tests/*.cpp 
	../tests/*.h)
add_executable (mast_tests ${mast_test_source_files})

target_link_libraries (mast_tests
			mast 
			${boost_unit_test_lib})
set_property (TARGET mast_tests APPEND 
	     PROPERTY INCLUDE_DIRECTORIES 
	     ${libmesh_dir}/include
	     ${}) 


ENABLE_TESTING()

ADD_TEST(NAME mast_tests
	 COMMAND mast_tests) 
	 #WORKING_DIRECTORY ../tests)

