#!/usr/bin/env bash

cd "${HOME}" || exit

if [ "${TRAVIS_OS_NAME}" = linux ]; then # Ubuntu Linux

  if [ "${TRAVIS_DIST}" = xenial ]; then # Ubuntu 16.04 Xenial Xerus
    wget -nv "https://github.com/MASTmultiphysics/mast-ci-packages/releases/download/libmesh-${LIBMESH_VERSION}-1.deb/libmesh-${LIBMESH_VERSION}-1.deb"
    sudo apt install "./libmesh-${LIBMESH_VERSION}-1.deb"

  # elif [ "${TRAVIS_DIST}" = bionic ]; then # Ubuntu 18.04 Bionic Beaver

  else
    echo "INVALID LINUX DISTRO: ${TRAVIS_DIST}"
    exit 1
  fi

elif [ "${TRAVIS_OS_NAME}" = osx ]; then # macOS 10.14, XCode 10.2
  mkdir "Code"
  mkdir "Code/spack-mstc"
  mkdir "Code/spack-mstc/spack"
  cd "Code/spack-mstc/spack" || exit
  wget -nv https://github.com/MASTmultiphysics/mast-ci-packages/releases/download/libmesh_macos_xcode12.2/libmesh_multiple_versions_macos_xcode_10.2.zip
  unzip -qq libmesh_multiple_versions_macos_xcode_10.2.zip

  #command -v opt/spack/darwin-mojave-x86_64/clang-10.0.1-apple/gcc-9.1.0-4amtftgtal2cnomzekpogzanzv6weadk/bin/gfortran
  #opt/spack/darwin-mojave-x86_64/clang-10.0.1-apple/gcc-9.1.0-4amtftgtal2cnomzekpogzanzv6weadk/bin/gfortran  --version
  #opt/spack/darwin-mojave-x86_64/clang-10.0.1-apple/gcc-9.1.0-4amtftgtal2cnomzekpogzanzv6weadk/bin/gfortran

else
  echo "INVALID OS: ${TRAVIS_OS_NAME}"
  exit 1
fi