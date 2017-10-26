MAST INSTALLATION
===============================

DEPENDENCY LIST
-------------------------------
MAST depends on the following libraries:

- PETSc (http://www.mcs.anl.gov/petsc/)
  MAST has been tested with PETSc version 3.6.3 and 3.7.7.

- SLEPc (http://slepc.upv.es)
  This builds on top of PETSc and provides the eigensolvers. The
  version numbers are in sync with that of PETSc, so a 3.6.x/3.7.x
  version of this library should be used.

- libMesh (http://libmesh.github.io)
  `git clone git://github.com/libMesh/libmesh.git`

- MPI

- LAPACK

- BLAS

- BOOST (http://www.boost.org)

- BOOST unit test framework

- EIGEN (http://eigen.tuxfamily.org/index.php?title=Main_Page)


- PARMETIS/METIS (http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview)
  Both PETSc and libMesh use these libraries. PETSc can be requested
  to download and install a version of this library during compile
  time. libMesh includes these libraries in the contrib
  subdirectory, which it can compile and link to.

- HDF5 (https://www.hdfgroup.org/HDF5/)
  This is used for the ExodusII output formats, which are efficiently
  read into Paraview. PETSc can be requested to download/build HDF5.
  libMesh also provides interfaces to Tecplot, and can be configured to
  build with Tecplot, in which case ExodusII becomes optional.

- GCMMA
  Presently, this is the default optimization library. Please obtain it
  from the author of GCMMA at: krille@math.kth.se

- libgfortran
  GCMMA is written in fortran. Hence, linking to GCMMA requires that
  MAST be linked to libgfortran.

CMAKE BUILD INSTRUCTIONS
-------------------------------

1. Download, build and install METIS and PARMETIS (PETSC can be request to
   install this during its configuration, which case skip this step.)
2. Download, build and install PETSc with MPI support. Please make
   sure to build the shared verion of this library.
3. Download, build and install SLEPc using the PETSc installation.
4. Download libMesh. Configure libMesh. Following configuration options
   are used locally. Please change the options to suit your local
   system. Note that shared library is being build.

```
PETSC_DIR=/Users/manav/Documents/codes/numerical_lib/petsc/\
SLEPC_DIR=/Users/manav/Documents/codes/numerical_lib/slepc/\
FC=mpif90-openmpi-mp F77=mpif90-openmpi-mp  CC=mpicc-openmpi-mp\
CXX=mpicxx-openmpi-mp ./configure\
--prefix=${PWD}/../ --enable-mpi --disable-unique-id\
--enable-dependency-tracking --enable-fortran --enable-shared\
--enable-exceptions --disable-openmp --disable-default-comm-world\
--enable-tracefiles  --enable-amr  --enable-vsmoother\
--enable-periodic  --enable-dirichlet  --enable-parmesh\
--enable-nodeconstraint  --enable-ghosted  --enable-pfem\
--enable-ifem --enable-second --enable-xdr --enable-reference-counting\
--enable-perflog --enable-examples --enable-boost --disable-trilinos\
--disable-tbb --enable-sfc --disable-tecplot --disable-tecio\
--enable-metis --enable-parmetis --enable-tetgen --enable-triangle\
--disable-vtk --enable-hdf5 --enable-libHilbert --enable-nanoflann\
--enable-exodus --enable-netcdf --enable-petsc --enable-slepc\
--with-mpi=/opt/local --with-metis=internal --with-hdf5=/opt/local/\
--with-methods="opt dbg"\
```

5. Build and install libMesh
6. Create a subdirectory MAST_DIR/build/opt to build MAST. In the
   subdirectory, configure cmake build system using the following
   command. Please modify it based on your local system
   configuration. You can use `-DCMAKE_BUILD_TYPE=Debug` to build the
   debug version.

```
cmake ../ -Dlibmesh_dir=~/Documents/codes/libmesh \
-Dpetsc_dir=~/Documents/codes/numerical_lib/petsc \
-Dslepc_dir=~/Documents/codes/numerical_lib/slepc \
-Dboost_include_dir=/opt/local/include \
-Dmpi_include_dir=/opt/local/include/openmpi-mp \
-Dmpi_lib_dir=/opt/local/lib/openmpi-mp -Dlapack_lib_dir=/usr/lib \
-Dblas_lib_dir=/usr/lib -Dboost_lib_dir=/opt/local/lib \
-Dgcmma_lib_file=gcmma \
-Dgcmma_lib_dir=~/Documents/codes/numerical_lib/gcmma/build \
-Ddot_lib_file=dot \
-Ddot_lib_dir=~/Documents/codes/numerical_lib/optimization_codes/dot/build \
-Dnpsol_lib_file=npsol \
-Dnpsol_lib_dir=~/Documents/codes/numerical_lib/optimization_codes/npsol/build \
-Dfortran_lib_file=libgfortran.3.dylib \
-Dfortran_lib_dir=/opt/local/lib/libgcc \
-Dboost_test_lib=boost_unit_test_framework-mt \
-DCMAKE_Fortran_COMPILER=mpif90 -DCMAKE_BUILD_TYPE=Release \
-DCMAKE_CXX_FLAGS=-std=c++11 \
-DENABLE_GCMMA=ON \
-DENABLE_DOT=ON \
-DENABLE_NPSOL=ON
```

Note: `CMAKE_BUILD_TYPE=Release` builds an optimized version, while
`CMAKE_BUILD_TYPE=Debug` builds the debug version of the code.

Note: GCMMA, DOT, NPSOL can be disabled by providing OFF as the compilation option
above. 

7. Build the library using
   `make mast`

8. Build one of the examples using.
   ` make example_driver`

9. Get instructions to run example by call it
   `./example_driver`


CLion IDE
-------------------------------
CLion is a C/C++ IDE produced by JetBrains (https://www.jetbrains.com/clion/)
that leverages the CMake build process for project organization.

To develop/build MAST using CLion, CMake options that are typically
supplied to the `cmake` terminal are provided in the preferences under
`Preferences > Build, Execution, Deployment > CMake` in the
`CMake Options:` box.

An example set of options is:

```
-DCMAKE_C_COMPILER=/usr/local/bin/mpicc
-DCMAKE_CXX_COMPILER=/usr/local/bin/mpicxx
-DCMAKE_FORTRAN_COMPILER=/usr/local/bin/mpifort
-Dlibmesh_dir=~/Code/libmesh-github-install
-Deigen_include_dir=/usr/local/include
-Dboost_include_dir=/usr/local/include
-Dboost_lib_dir=/usr/local/lib
-Dboost_test_lib=boost_unit_test_framework
-Dmpi_include_dir=/usr/local/include
-Dmpi_lib_dir=/usr/local/lib
-Dblas_lib_dir=/usr/lib
-Dlapack_lib_dir=/usr/lib
-Dpetsc_dir=~/Code/mast-multiphysics-deps2
-Dslepc_dir=~/Code/mast-multiphysics-deps2
-Dfortran_lib_dir=/usr/local/lib/gcc/7
-Dfortran_lib_file=libgfortran.dylib
```

You will want to use options corresponding to your own environment as
described in the *CMAKE BUILD INSTRUCTIONS* section.