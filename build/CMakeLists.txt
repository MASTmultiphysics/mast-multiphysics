project (MAST)

cmake_minimum_required(VERSION 2.8)

enable_language (CXX Fortran)

set (MAST_VERSION_MAJOR 0)
set (MAST_VERSION_MINOR 9)

set (libmesh_dir       "libmesh-dir" CACHE STRING "Directory containing libMesh header and library")
set (vtk_include_dir   "vtk-include-dir" CACHE STRING "Directory containing VTK header files")
set (vtk_lib_dir       "vtk-lib-dir" CACHE STRING "Directory containing  VTK library files" )
set (vtk_lib_file      "vtkIO vtkCommon vtkFiltering" CACHE STRING "VTK library files" )
separate_arguments (vtk_lib_file)
set (tbb_include_dir   "tbb-include-dir" CACHE STRING "Directory containing TBB header files")
set (tbb_lib_dir       "tbb-lib-dir" CACHE STRING "Directory containing TBB library files")
set (tbb_lib_file      "tbb tbbmalloc" CACHE STRING "TBB library files")
separate_arguments (tbb_lib_file)
set (boost_include_dir "boost-include-dir" CACHE STRING "Directory containing boost header files")
set (boost_test_lib    "boost-unit-test-lib" CACHE STRING "BOOST unit test library")
set (boost_lib_dir     "boost-lib-dir" CACHE STRING "Directory containing boost unit test library files")
set (mpi_include_dir   "mpi-include-dir" CACHE STRING "Directory containing MPI header files")
set (mpi_lib_dir       "mpi-lib-dir" CACHE STRING "Directory containing MPI library files")
set (lapack_lib_file   "lapack" CACHE STRING "LAPACK library")
separate_arguments (lapack_lib_file)
set (blas_lib_file     "blas" CACHE STRING "BLAS library")
separate_arguments (blas_lib_file)
set (lapack_lib_dir    "lapack-lib-dir" CACHE STRING "Directory containing LAPACK library")
set (blas_lib_dir      "blas-lib-dir" CACHE STRING "Directory containing BLAS library")
set (CMAKE_Fortran_COMPILER "fortran-compiler" CACHE STRING "Fortran compiler")


#set (EXTRA_LIBS ${EXTRA_LIBS} MathFunctions)


include_directories (../src 
                     ../examples
		    ${vtk_include_dir}
		    ${tbb_include_dir}
		    ${boost_include_dir})

include_directories (SYSTEM ${mpi_include_dir})


find_library (mesh_dbg mesh_dbg ${libmesh_dir}/lib)
find_library (mesh_opt mesh_opt ${libmesh_dir}/lib)

find_library (mpicxx_lib    mpi_cxx              ${mpi_lib_dir})
find_library (mpi_lib       mpi                  ${mpi_lib_dir})

find_library (boost_unit_test_lib   ${boost_test_lib}    ${boost_lib_dir})

foreach (f ${tbb_lib_file})
   find_library (tbb_lib    ${f}   ${tbb_lib_dir})
endforeach (f)

foreach (f ${vtk_lib_file})
   find_library (vtk_lib    ${f}   ${vtk_lib_dir})
endforeach (f)

foreach (f ${lapack_lib_file})
   find_library (lapack_lib    ${f}   ${lapack_lib_dir})
endforeach (f)

foreach (f ${blas_lib_file})
find_library (blas_lib      ${f}     ${blas_lib_dir})
endforeach (f)


#
#  add files for GCMMA library and add it to cmake's build list
#
#file (GLOB gcmma_source_files ../src/Optimization/gcmma/*.f)
#add_library (gcmma ${gcmma_source_files})


#
#  add files for MAST library and add it to cmake's build list
#
file (GLOB_RECURSE mast_source_files ../src/*.cpp ../src/*.h)
add_library (mast ${mast_source_files})

target_link_libraries (mast  
		      debug ${mesh_dbg}
		      optimized ${mesh_opt}
		      #gcmma
		      ${lapack_lib}
		      ${blas_lib}
		      ${vtk_lib}
		      ${tbb_lib}
		      ${mpicxx_lib}
		      ${mpi_lib})
set_property (TARGET mast APPEND 
	     PROPERTY INCLUDE_DIRECTORIES 
	     ${libmesh_dir}/include) 


#
#  tell cmake to link gaussian_bump_convergence example with mast library
#
file (GLOB_RECURSE ex_1_source_files 
	../examples/ex_1/*.cpp 
	../examples/ex_1/*.h)
add_executable (ex_1 ${ex_1_source_files})

target_link_libraries (ex_1
			mast)
set_property (TARGET ex_1 APPEND 
	     PROPERTY INCLUDE_DIRECTORIES 
	     ${libmesh_dir}/include) 


