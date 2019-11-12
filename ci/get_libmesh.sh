#!/usr/bin/env bash

cd ${HOME}

if [ "${TRAVIS_OS_NAME}" = linux ]; then # Ubuntu Linux

  if [ "${TRAVIS_DIST}" = xenial ]; then # Ubuntu 16.04 Xenial Xerus
    wget -nv https://github.com/MASTmultiphysics/mast-ci-packages/releases/download/libmesh-${LIBMESH_VERSION}-1.deb/libmesh-${LIBMESH_VERSION}-1.deb
    sudo apt install ./libmesh-${LIBMESH_VERSION}-1.deb

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

elif [ "${TRAVIS_OS_NAME}" = osx ]; then # macOS 10.14, XCode
  echo "Hello From OSX"

else
  echo "INVALID OS: ${TRAVIS_OS_NAME}"
  exit 1
fi