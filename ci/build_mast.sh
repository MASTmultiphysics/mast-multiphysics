#!/usr/bin/env bash

# Steps common to all OS/toolchains.
export MAST_INSTALL_DIR=${HOME}/mast || exit
cd "${TRAVIS_BUILD_DIR}" || exit

# Create directories to test both Release (optimized) and Debug build.
mkdir build_rel || exit
mkdir build_dbg || exit
export REL_BUILD_DIR=${TRAVIS_BUILD_DIR}/build_rel
export DBG_BUILD_DIR=${TRAVIS_BUILD_DIR}/build_dbg

# Outer if chooses Linux/macOS based on Travis CI worker.
if [ "${TRAVIS_OS_NAME}" = linux ]; then # Ubuntu Linux

  if [ "${TRAVIS_DIST}" = xenial ]; then # Ubuntu 16.04 Xenial Xerus

    # First let us build/install a Debug version of MAST (-DCMAKE_BUILD_TYPE=Debug).
    echo "TEST DEBUG BUILD..."

    if [ "${LIBMESH_VERSION}" = "1.3.1" ]; then # No Debug build for libMesh < 1.5.0
      echo "-- NO DEBUG BUILD WITH libMesh < 1.5.0"
    elif [ "${LIBMESH_VERSION}" = "1.4.1" ]; then
      echo "-- NO DEBUG BUILD WITH libMesh < 1.5.0"
    else

      cd "${DBG_BUILD_DIR}" || exit
      /opt/cmake/bin/cmake .. \
        -DCMAKE_BUILD_TYPE=Debug \
        -DCMAKE_INSTALL_PREFIX="${MAST_INSTALL_DIR}" \
        -DCMAKE_C_COMPILER=mpicc \
        -DCMAKE_CXX_COMPILER=mpic++ \
        -DCMAKE_Fortran_COMPILER=mpifort \
        -DlibMesh_DIR=/usr/local \
        -DPETSc_DIR=/usr/lib/petscdir/3.6.2/x86_64-linux-gnu-real \
        -DSLEPc_DIR=/usr/lib/slepcdir/3.6.1/x86_64-linux-gnu-real \
        -DEIGEN3_ROOT=/usr/include/eigen3 \
        -DPython3_DIR=/usr \
        -DBOOST_ROOT=/usr \
        -DBUILD_DOC=ON \
        -DENABLE_DOT=OFF \
        -DENABLE_GCMMA=OFF \
        -DENABLE_SNOPT=OFF\
        -DENABLE_NASTRANIO=ON \
        -DENABLE_CYTHON=ON || exit

      if [ ${CI_BUILD_DOCS} ]; then
        echo "No CI documentation for a Debug build."
      else
        make -j 2 || exit
        echo "RUNNING UNIT TESTS" || exit
        cd "${DBG_BUILD_DIR}/tests" || exit
        ctest --force-new-ctest-process --output-on-failure --timeout 10
        echo "RUNNING SHORT EXAMPLES" || exit
        cd "${DBG_BUILD_DIR}/examples"  || exit
        ctest --force-new-ctest-process --output-on-failure -L "SHORT" --timeout 60
      fi

    fi # No Debug build for libMesh < 1.5.0

    # Now build/install a Release (optimized) version of MAST (-DCMAKE_BUILD_TYPE=Release).
    echo "TEST RELEASE/OPTIMIZED BUILD..."
    cd "${REL_BUILD_DIR}" || exit
    /opt/cmake/bin/cmake .. \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX="${MAST_INSTALL_DIR}" \
      -DCMAKE_C_COMPILER=mpicc \
      -DCMAKE_CXX_COMPILER=mpic++ \
      -DCMAKE_Fortran_COMPILER=mpifort \
      -DlibMesh_DIR=/usr/local \
      -DPETSc_DIR=/usr/lib/petscdir/3.6.2/x86_64-linux-gnu-real \
      -DSLEPc_DIR=/usr/lib/slepcdir/3.6.1/x86_64-linux-gnu-real \
      -DEIGEN3_ROOT=/usr/include/eigen3 \
      -DPython3_DIR=/usr \
      -DBOOST_ROOT=/usr \
      -DBUILD_DOC=ON \
      -DENABLE_DOT=OFF \
      -DENABLE_GCMMA=OFF \
      -DENABLE_SNOPT=OFF \
      -DENABLE_NASTRANIO=ON \
      -DENABLE_CYTHON=ON || exit

    if [ ${CI_BUILD_DOCS} ]; then
      make doc_doxygen || exit
      cd "${TRAVIS_BUILD_DIR}" || exit
      ci/prepare_docs.sh || exit
    else
      make -j 2 || exit
      make install || exit
      echo "RUNNING UNIT TESTS" || exit
      cd "${REL_BUILD_DIR}/tests" || exit
      ctest --force-new-ctest-process --output-on-failure --timeout 10
      echo "RUNNING SHORT EXAMPLES" || exit
      cd "${REL_BUILD_DIR}/examples"  || exit
      ctest --force-new-ctest-process --output-on-failure -L "SHORT" --timeout 60
    fi

  elif [ "${TRAVIS_DIST}" = bionic ]; then # Ubuntu 18.04 Bionic Beaver

    # First let us build/install a Debug version of MAST (-DCMAKE_BUILD_TYPE=Debug).
    echo "TEST DEBUG BUILD..."

    if [ "${LIBMESH_VERSION}" = "1.3.1" ]; then # No Debug build for libMesh < 1.5.0
      echo "-- NO DEBUG BUILD WITH libMesh < 1.5.0"
    elif [ "${LIBMESH_VERSION}" = "1.4.1" ]; then
      echo "-- NO DEBUG BUILD WITH libMesh < 1.5.0"
    else

      cd "${DBG_BUILD_DIR}" || exit
      /opt/cmake/bin/cmake .. \
        -DCMAKE_BUILD_TYPE=Debug \
        -DCMAKE_INSTALL_PREFIX="${MAST_INSTALL_DIR}" \
        -DCMAKE_C_COMPILER=mpicc \
        -DCMAKE_CXX_COMPILER=mpic++ \
        -DCMAKE_Fortran_COMPILER=mpifort \
        -DlibMesh_DIR=/usr/local \
        -DPETSc_DIR=/usr/lib/petsc \
        -DSLEPc_DIR=/usr/lib/slepc \
        -DEIGEN3_ROOT=/usr/include/eigen3 \
        -DPython3_DIR=/usr \
        -DBOOST_ROOT=/usr \
        -DBUILD_DOC=ON \
        -DENABLE_DOT=OFF \
        -DENABLE_GCMMA=OFF \
        -DENABLE_SNOPT=OFF\
        -DENABLE_NASTRANIO=ON \
        -DENABLE_CYTHON=ON || exit

      if [ ${CI_BUILD_DOCS} ]; then
        echo "No CI documentation for a Debug build."
      else
        make -j 2 || exit
        echo "RUNNING UNIT TESTS" || exit
        cd "${DBG_BUILD_DIR}/tests" || exit
        ctest --force-new-ctest-process --output-on-failure --timeout 10
        echo "RUNNING SHORT EXAMPLES" || exit
        cd "${DBG_BUILD_DIR}/examples"  || exit
        ctest --force-new-ctest-process --output-on-failure -L "SHORT" --timeout 60
      fi

    fi # No Debug build for libMesh < 1.5.0

    # Now build/install a Release (optimized) version of MAST (-DCMAKE_BUILD_TYPE=Release).
    echo "TEST RELEASE/OPTIMIZED BUILD..."
    cd "${REL_BUILD_DIR}" || exit
    /opt/cmake/bin/cmake .. \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX="${MAST_INSTALL_DIR}" \
      -DCMAKE_C_COMPILER=mpicc \
      -DCMAKE_CXX_COMPILER=mpic++ \
      -DCMAKE_Fortran_COMPILER=mpifort \
      -DlibMesh_DIR=/usr/local \
      -DPETSc_DIR=/usr/lib/petsc \
      -DSLEPc_DIR=/usr/lib/slepc \
      -DEIGEN3_ROOT=/usr/include/eigen3 \
      -DPython3_DIR=/usr \
      -DBOOST_ROOT=/usr \
      -DBUILD_DOC=ON \
      -DENABLE_DOT=OFF \
      -DENABLE_GCMMA=OFF \
      -DENABLE_SNOPT=OFF \
      -DENABLE_NASTRANIO=ON \
      -DENABLE_CYTHON=ON || exit

    if [ ${CI_BUILD_DOCS} ]; then
      make doc_doxygen || exit
      cd "${TRAVIS_BUILD_DIR}" || exit
      ci/prepare_docs.sh || exit
    else
      make -j 2 || exit
      make install || exit
      echo "RUNNING UNIT TESTS" || exit
      cd "${REL_BUILD_DIR}/tests" || exit
      ctest --force-new-ctest-process --output-on-failure --timeout 10
      echo "RUNNING SHORT EXAMPLES" || exit
      cd "${REL_BUILD_DIR}/examples"  || exit
      ctest --force-new-ctest-process --output-on-failure -L "SHORT" --timeout 60
    fi


  else
    echo "INVALID LINUX DISTRO: ${TRAVIS_DIST}"
    exit 1
  fi

elif [ "${TRAVIS_OS_NAME}" = osx ]; then # macOS 10.14, XCode 10.2

  if [ "${LIBMESH_VERSION}" = "1.3.1" ]; then
    PETSc_DIR=/Users/travis/Code/spack-mstc/spack/opt/spack/darwin-mojave-x86_64/clang-10.0.1-apple/petsc-3.11.4-4ode2ljnkurc55362bc6k6wtlgfjpdmf
    SLEPc_DIR=/Users/travis/Code/spack-mstc/spack/opt/spack/darwin-mojave-x86_64/clang-10.0.1-apple/slepc-3.11.2-nts2wrpmixfsvfev7x5b5e2e4vpct6wy
    libMesh_DIR=/Users/travis/Code/spack-mstc/spack/opt/spack/darwin-mojave-x86_64/clang-10.0.1-apple/libmesh-1.3.1-4h7zpmhjpw6qwodrpoehwaq7fukjwfno
  elif [ "${LIBMESH_VERSION}" = "1.4.1" ]; then
    PETSc_DIR=/Users/travis/Code/spack-mstc/spack/opt/spack/darwin-mojave-x86_64/clang-10.0.1-apple/petsc-3.11.4-4ode2ljnkurc55362bc6k6wtlgfjpdmf
    SLEPc_DIR=/Users/travis/Code/spack-mstc/spack/opt/spack/darwin-mojave-x86_64/clang-10.0.1-apple/slepc-3.11.2-nts2wrpmixfsvfev7x5b5e2e4vpct6wy
    libMesh_DIR=/Users/travis/Code/spack-mstc/spack/opt/spack/darwin-mojave-x86_64/clang-10.0.1-apple/libmesh-1.4.1-bry6shpqjekanfoepltqyic5ami5umor
  elif [ "${LIBMESH_VERSION}" = "1.5.0" ]; then
    PETSc_DIR=/Users/travis/Code/spack-mstc/spack/opt/spack/darwin-mojave-x86_64/clang-10.0.1-apple/petsc-3.12.1-edmddchkkhtwckfumww3mtghyadytdxj
    SLEPc_DIR=/Users/travis/Code/spack-mstc/spack/opt/spack/darwin-mojave-x86_64/clang-10.0.1-apple/slepc-3.12.0-dbsf3rtgj4hoih6okdja6godorotij22
    libMesh_DIR=/Users/travis/Code/spack-mstc/spack/opt/spack/darwin-mojave-x86_64/clang-10.0.1-apple/libmesh-1.5.0-dcgl6zldkp4xvfvfrajghwoy5b24ilrr
  elif [ "${LIBMESH_VERSION}" = "1.5.1" ]; then
    PETSc_DIR=/Users/travis/Code/spack-mstc/spack/opt/spack/darwin-mojave-x86_64/clang-10.0.1-apple/petsc-3.12.1-edmddchkkhtwckfumww3mtghyadytdxj
    SLEPc_DIR=/Users/travis/Code/spack-mstc/spack/opt/spack/darwin-mojave-x86_64/clang-10.0.1-apple/slepc-3.12.0-dbsf3rtgj4hoih6okdja6godorotij22
    libMesh_DIR=/Users/travis/Code/spack-mstc/spack/opt/spack/darwin-mojave-x86_64/clang-10.0.1-apple/libmesh-1.5.1-5y4kmrrnq6axqa2hv2ag36tjv6mgwknw
  fi

  # First let us build/install a Debug version of MAST (-DCMAKE_BUILD_TYPE=Debug).
  echo "TEST DEBUG BUILD..."
  cd "${DBG_BUILD_DIR}" || exit
  /Users/travis/Code/spack-mstc/spack/opt/spack/darwin-mojave-x86_64/clang-10.0.1-apple/cmake-3.14.4-vlzpkvowre6hweutogn563ouzkb73jsq/bin/cmake .. \
    -DCMAKE_BUILD_TYPE=Debug \
    -DCMAKE_INSTALL_PREFIX="${MAST_INSTALL_DIR}" \
    -DCMAKE_C_COMPILER=/Users/travis/Code/spack-mstc/spack/opt/spack/darwin-mojave-x86_64/clang-10.0.1-apple/openmpi-3.1.4-sep4omvhokkexyojwrahdckfugactovb/bin/mpicc \
    -DCMAKE_CXX_COMPILER=/Users/travis/Code/spack-mstc/spack/opt/spack/darwin-mojave-x86_64/clang-10.0.1-apple/openmpi-3.1.4-sep4omvhokkexyojwrahdckfugactovb/bin/mpic++ \
    -DCMAKE_Fortran_COMPILER=/Users/travis/Code/spack-mstc/spack/opt/spack/darwin-mojave-x86_64/clang-10.0.1-apple/openmpi-3.1.4-sep4omvhokkexyojwrahdckfugactovb/bin/mpif90 \
    -DMPIEXEC_EXECUTABLE=/Users/travis/Code/spack-mstc/spack/opt/spack/darwin-mojave-x86_64/clang-10.0.1-apple/openmpi-3.1.4-sep4omvhokkexyojwrahdckfugactovb/bin/mpiexec \
    -DlibMesh_DIR="${libMesh_DIR}" \
    -DPETSc_DIR="${PETSc_DIR}" \
    -DSLEPc_DIR="${SLEPc_DIR}" \
    -DBOOST_ROOT=/Users/travis/Code/spack-mstc/spack/opt/spack/darwin-mojave-x86_64/clang-10.0.1-apple/boost-1.70.0-3mel532oks2cxz76mqrxvn4xrl3bg5ri \
    -DHDF5_ROOT=/Users/travis/Code/spack-mstc/spack/opt/spack/darwin-mojave-x86_64/clang-10.0.1-apple/hdf5-1.10.5-kfyeqhtsenqumr7eut6b47mjev7tstom \
    -DEIGEN3_ROOT=/Users/travis/Code/spack-mstc/spack/opt/spack/darwin-mojave-x86_64/clang-10.0.1-apple/eigen-3.3.7-bd2r4aqrkox7dpebj2r3gqvgpqzwuh7x \
    -DLAPACK_LIBRARIES=/Users/travis/Code/spack-mstc/spack/opt/spack/darwin-mojave-x86_64/clang-10.0.1-apple/openblas-0.3.7-xqeap7iegoomce3es67cd7exlnq3neue/lib/libopenblas.dylib \
    -DBLAS_LIBRARIES=/Users/travis/Code/spack-mstc/spack/opt/spack/darwin-mojave-x86_64/clang-10.0.1-apple/openblas-0.3.7-xqeap7iegoomce3es67cd7exlnq3neue/lib/libopenblas.dylib \
    -DPython3_DIR=/Users/travis/Code/spack-mstc/spack/opt/spack/darwin-mojave-x86_64/clang-10.0.1-apple/python-3.7.4-hh7odiqwzpw2nnfdcfzjxvb3ggyzwvfk \
    -DENABLE_GCMMA=OFF \
    -DENABLE_DOT=OFF \
    -DENABLE_SNOPT=OFF \
    -DENABLE_NLOPT=OFF \
    -DENABLE_CYTHON=OFF \
    -DBUILD_DOC=OFF \
    -DENABLE_NASTRANIO=ON \
    -DENABLE_CYTHON=ON || exit

  make -j 2 || exit
  echo "RUNNING UNIT TESTS" || exit
  cd "${DBG_BUILD_DIR}/tests" || exit
  ctest --force-new-ctest-process --output-on-failure --timeout 10
  echo "RUNNING SHORT EXAMPLES" || exit
  cd "${DBG_BUILD_DIR}/examples" || exit
  ctest --force-new-ctest-process --output-on-failure -L "SHORT" --timeout 60

  # Now build/install a Release (optimized) version of MAST (-DCMAKE_BUILD_TYPE=Release).
  echo "TEST RELEASE/OPTIMIZED BUILD..."
  cd "${REL_BUILD_DIR}" || exit
  /Users/travis/Code/spack-mstc/spack/opt/spack/darwin-mojave-x86_64/clang-10.0.1-apple/cmake-3.14.4-vlzpkvowre6hweutogn563ouzkb73jsq/bin/cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="${MAST_INSTALL_DIR}" \
    -DCMAKE_C_COMPILER=/Users/travis/Code/spack-mstc/spack/opt/spack/darwin-mojave-x86_64/clang-10.0.1-apple/openmpi-3.1.4-sep4omvhokkexyojwrahdckfugactovb/bin/mpicc \
    -DCMAKE_CXX_COMPILER=/Users/travis/Code/spack-mstc/spack/opt/spack/darwin-mojave-x86_64/clang-10.0.1-apple/openmpi-3.1.4-sep4omvhokkexyojwrahdckfugactovb/bin/mpic++ \
    -DCMAKE_Fortran_COMPILER=/Users/travis/Code/spack-mstc/spack/opt/spack/darwin-mojave-x86_64/clang-10.0.1-apple/openmpi-3.1.4-sep4omvhokkexyojwrahdckfugactovb/bin/mpif90 \
    -DMPIEXEC_EXECUTABLE=/Users/travis/Code/spack-mstc/spack/opt/spack/darwin-mojave-x86_64/clang-10.0.1-apple/openmpi-3.1.4-sep4omvhokkexyojwrahdckfugactovb/bin/mpiexec \
    -DlibMesh_DIR="${libMesh_DIR}" \
    -DPETSc_DIR="${PETSc_DIR}" \
    -DSLEPc_DIR="${SLEPc_DIR}" \
    -DBOOST_ROOT=/Users/travis/Code/spack-mstc/spack/opt/spack/darwin-mojave-x86_64/clang-10.0.1-apple/boost-1.70.0-3mel532oks2cxz76mqrxvn4xrl3bg5ri \
    -DHDF5_ROOT=/Users/travis/Code/spack-mstc/spack/opt/spack/darwin-mojave-x86_64/clang-10.0.1-apple/hdf5-1.10.5-kfyeqhtsenqumr7eut6b47mjev7tstom \
    -DEIGEN3_ROOT=/Users/travis/Code/spack-mstc/spack/opt/spack/darwin-mojave-x86_64/clang-10.0.1-apple/eigen-3.3.7-bd2r4aqrkox7dpebj2r3gqvgpqzwuh7x \
    -DLAPACK_LIBRARIES=/Users/travis/Code/spack-mstc/spack/opt/spack/darwin-mojave-x86_64/clang-10.0.1-apple/openblas-0.3.7-xqeap7iegoomce3es67cd7exlnq3neue/lib/libopenblas.dylib \
    -DBLAS_LIBRARIES=/Users/travis/Code/spack-mstc/spack/opt/spack/darwin-mojave-x86_64/clang-10.0.1-apple/openblas-0.3.7-xqeap7iegoomce3es67cd7exlnq3neue/lib/libopenblas.dylib \
    -DPython3_DIR=/Users/travis/Code/spack-mstc/spack/opt/spack/darwin-mojave-x86_64/clang-10.0.1-apple/python-3.7.4-hh7odiqwzpw2nnfdcfzjxvb3ggyzwvfk \
    -DENABLE_GCMMA=OFF \
    -DENABLE_DOT=OFF \
    -DENABLE_SNOPT=OFF \
    -DENABLE_NLOPT=OFF \
    -DENABLE_CYTHON=OFF \
    -DBUILD_DOC=OFF \
    -DENABLE_NASTRANIO=ON \
    -DENABLE_CYTHON=ON || exit

  make -j 2 || exit
  make install || exit

  echo "RUNNING UNIT TESTS" || exit
  cd "${REL_BUILD_DIR}/tests" || exit
  ctest --force-new-ctest-process --output-on-failure --timeout 10
  echo "RUNNING SHORT EXAMPLES" || exit
  cd "${REL_BUILD_DIR}/examples" || exit
  ctest --force-new-ctest-process --output-on-failure -L "SHORT" --timeout 60

else
  echo "INVALID OS: ${TRAVIS_OS_NAME}"
  exit 1
fi

cd "${TRAVIS_BUILD_DIR}" || exit