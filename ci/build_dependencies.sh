#!/usr/bin/env bash

if [ "${MY_OS_TYPE}" = Ubuntu ]; then # Ubuntu Linux

  if [ "${MY_OS}" = ubuntu-20.04 ]; then # Ubuntu 20.04 Focal Fossa
    # Regular libMesh/MAST dependencies.
    sudo apt-get update
    sudo apt-get -qq install -y \
      build-essential gfortran wget less m4 git cmake \
      python3-all python3-all-dev python3-all-dbg \
      openmpi-bin libopenmpi-dev python3-mpi4py python3-mpi4py-dbg \
      petsc-dev python3-petsc4py \
      slepc-dev python3-slepc4py \
      libparpack2-dev \
      libmetis-dev \
      libnetcdf-dev \
      libboost-all-dev \
      libeigen3-dev \
      libnlopt-dev \
      libadolc-dev \
      doxygen graphviz rsync

  elif [ "${MY_OS}" = ubuntu-18.04 ]; then # Ubuntu 18.04 Bionic Beaver
    # Regular libMesh/MAST dependencies.
    sudo apt-get update
      sudo apt-get -qq install -y \
        build-essential gfortran wget less m4 git cmake \
        python3-all python3-all-dev python3-all-dbg \
        openmpi-bin libopenmpi-dev python3-mpi4py \
        petsc-dev python3-petsc4py \
        slepc-dev python3-slepc4py \
        libparpack2-dev \
        libmetis-dev \
        libnetcdf-dev \
        libboost-all-dev \
        libeigen3-dev \
        libnlopt-dev \
        libadolc-dev \
        doxygen graphviz rsync
        # sudo apt-get -qq install -y texlive-latex-base dvi2ps ghostscript

  else
    echo "INVALID LINUX DISTRO: ${MY_OS}"
    exit 1
  fi

  # Get pip working with external Python 3.7.
  # wget https://bootstrap.pypa.io/get-pip.py || exit
  # sudo python3.7 get-pip.py || exit

  # sudo python3.7 -m pip install numpy scipy docopt colorama pandas h5py matplotlib cpylog pyNastran
  # sudo python3.7 -m pip install Cython --install-option="--no-cython-compile"

  # # Update to later CMake release.
  # wget https://github.com/Kitware/CMake/releases/download/v3.15.5/cmake-3.15.5-Linux-x86_64.sh || exit
  # sudo mkdir /opt/cmake || exit
  # sudo sh cmake-3.15.5-Linux-x86_64.sh --prefix=/opt/cmake --skip-license || exit

# elif [ "${TRAVIS_OS_NAME}" = osx ]; then # macOS 10.14, XCode 10.2
#   # Currently we don't do anything here since we get all dependencies for macOS
#   # from the binary download with "ci/get_libmesh.sh" in the next stage.
#   echo "Hello From OSX"

else
  echo "INVALID OS (MY_OS_TYPE): ${MY_OS_TYPE}"
  exit 1
fi