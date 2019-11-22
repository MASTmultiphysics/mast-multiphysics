#!/usr/bin/env bash

if [ "${TRAVIS_OS_NAME}" = linux ]; then # Ubuntu Linux

  if [ "${TRAVIS_DIST}" = xenial ]; then # Ubuntu 16.04 Xenial Xerus
    sudo apt-get -qq update
    sudo apt-get -qq install -y gfortran wget m4
    sudo apt-get -qq install -y openmpi-bin libopenmpi-dev
    sudo apt-get -qq install -y libpetsc3.6 libpetsc3.6.2-dev
    sudo apt-get -qq install -y libslepc3.6 libslepc3.6.1-dev libparpack2-dev
    sudo apt-get -qq install -y libboost-all-dev
    sudo apt-get -qq install -y libeigen3-dev
    sudo apt-get -qq install -y doxygen graphviz rsync
    sudo apt-get -qq install -y texlive-latex-base dvi2ps ghostscript

    # Update to later CMake release.
    cd ${HOME} || exit
    wget https://github.com/Kitware/CMake/releases/download/v3.15.5/cmake-3.15.5-Linux-x86_64.sh || exit
    sudo mkdir /opt/cmake || exit
    sudo sh cmake-3.15.5-Linux-x86_64.sh --prefix=/opt/cmake || exit

  # elif [ "${TRAVIS_DIST}" = bionic ]; then # Ubuntu 18.04 Bionic Beaver
  #  sudo apt-get -qq update
  #  sudo apt-get -qq install -y gfortran wget m4
  #  which gfortran
  #  gcc --version
  #  echo "Hello From BIONIC"

  else
    echo "INVALID LINUX DISTRO: ${TRAVIS_DIST}"
    exit 1
  fi

elif [ "${TRAVIS_OS_NAME}" = osx ]; then # macOS 10.14, XCode 10.2
  # Currently we don't do anything here since we get all dependencies for macOS
  # from the binary download with "ci/get_libmesh.sh" in the next stage.
  echo "Hello From OSX"

else
  echo "INVALID OS: ${TRAVIS_OS_NAME}"
  exit 1
fi