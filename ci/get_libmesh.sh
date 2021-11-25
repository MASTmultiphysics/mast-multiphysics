#!/usr/bin/env bash

cd "${MY_WORK_PATH}" || exit

if [ "${MY_OS_TYPE}" = Ubuntu ]; then # Ubuntu Linux

  if [ "${MY_OS}" = ubuntu-20.04 ]; then # Ubuntu 20.04 Focal Fossa
    wget -nv "https://github.com/MASTmultiphysics/mast-ci-packages/releases/download/gha-ubuntu20.04/libmesh-${MY_LIBMESH_VERSION}-1.deb" || exit
    sudo apt install "./libmesh-${MY_LIBMESH_VERSION}-1.deb" || exit

  elif [ "${MY_OS}" = ubuntu-18.04 ]; then # Ubuntu 18.04 Bionic Beaver
    wget -nv "https://github.com/MASTmultiphysics/mast-ci-packages/releases/download/gha-ubuntu18.04/libmesh-${MY_LIBMESH_VERSION}-1.deb" || exit
    sudo apt install "./libmesh-${MY_LIBMESH_VERSION}-1.deb" || exit

  else
    echo "INVALID LINUX DISTRO: ${MY_OS}"
    exit 1
  fi

  sudo ldconfig # Don't remove this, it is necessary to update the cache of shared libraries to those that were just
                # unpacked in the previous steps. For some reason `apt install` doesn't do this.

# elif [ "${TRAVIS_OS_NAME}" = osx ]; then # macOS 10.14, XCode 10.2
#   mkdir "Code" || exit
#   mkdir "Code/spack-mstc" || exit
#   mkdir "Code/spack-mstc/spack" || exit
#   cd "Code/spack-mstc/spack" || exit
#   wget -nv https://github.com/MASTmultiphysics/mast-ci-packages/releases/download/libmesh_macos_xcode12.2/libmesh_multiple_versions_macos_xcode_10.2.zip || exit
#   unzip -qq libmesh_multiple_versions_macos_xcode_10.2.zip || exit

  #command -v opt/spack/darwin-mojave-x86_64/clang-10.0.1-apple/gcc-9.1.0-4amtftgtal2cnomzekpogzanzv6weadk/bin/gfortran
  #opt/spack/darwin-mojave-x86_64/clang-10.0.1-apple/gcc-9.1.0-4amtftgtal2cnomzekpogzanzv6weadk/bin/gfortran  --version
  #opt/spack/darwin-mojave-x86_64/clang-10.0.1-apple/gcc-9.1.0-4amtftgtal2cnomzekpogzanzv6weadk/bin/gfortran

else
  echo "INVALID OS: ${MY_OS_TYPE}"
  exit 1
fi