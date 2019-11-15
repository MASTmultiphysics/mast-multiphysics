# Travis CI Testing/Deployment

On the GitHub repository of MAST, Travis CI (https://travis-ci.com) is utilized for build testing on multiple OS's and
deployment of updated documentation online. In the future we will integrate the execution of both unit tests as well as 
example problems for integration testing.

This directory contains scripts utilized by the Travis CI process.

## Travis CI Processes
When a commit is pushed to the MASTmultiphysics/mast-multiphysics GitHub repository, GitHub initiates the continuous
integration processes on Travis CI. What these processes do is controlled by the `.travis.yml` file in the root of the
repository. The different phases of a Travis CI are described at https://docs.travis-ci.com/user/job-lifecycle/. This
repository utilizes a relatively basic subset of the phases and overall capability of Travis CI.

The files in this directory are utilized in the main phases of the Travis CI lifecycle:

1. `before_install` - `build_dependencies.sh` - This script ensures that the Travis CI runner executing the current job
has the required dependencies for libMesh and building MAST.
2. `install` - `get_libmesh.sh` - This script fetches pre-built binaries of libMesh on Linux or the entire dependency
package on macOS. Travis CI jobs have a limited time for execution that is not long enough to build libMesh 
(which is not available in job runners standard Linux/macOS package managers.)
3. `script` - `build_mast.sh` - Depending on the job, this script either builds the MAST library/examples or builds the
doxygen documentation.
4. `deploy` - `deploy_docs.sh` - For the job satisfying the appropriate conditions, this script pushes the documentation
website produced by doxygen to the hosting location.

## Multiple OS's/Jobs/Environments
The `matrix` section of the `.travis.yml` file allows for the definition of multiple jobs to be run that can be done
in different environments or on different workers. Each item under `matrix/include` causes the execution of different 
Travis CI worker with the specified environment. In the current setup, we utilize both Linux and macOS to test builds 
against different compilers as well as different versions of the libMesh dependency, which is periodically updated by 
that projects maintainers. In this setup, Travis CI allows for the specification of system environment variables that
are utilized by the scripts in the `ci` folder for distinguishing work that is specific to each environment.

### Linux
The current Linux build environment utilizes Ubuntu 16.04. Dependencies available in the Ubuntu apt package repositories
are leveraged; however, in the future these may be updated as versions of PETSc/SLEPc are quite old. The current
compiler on Linux is GNU GCC-5.4.0. The libMesh dependencies are provided as external binaries that were built in an
identical environment and archived.

### macOS
The macOS build on Travis CI current utilizes macOS version 10.14.4 with the default clang C/C++ compiler (Apple LLVM 
version 10.0.1). Travis CI suggests utilizing Homebrew to install dependencies, which currently does not contain many of
the packages required for libMesh/MAST. To overcome this, a complete set dependencies (including multiple versions of
libMesh) is built/archived using Spack on the same environment setup (macOS version/compiler) the Travis CI utilizes. 
These dependencies are fetched by the Travis CI runner and MAST is built against them. Since macOS does not provide a
Fortran compiler, GCC-9.1.0 is first built via Spack. The compiler toolchain is then macOS-provided C/C++ and
gfortran-9.