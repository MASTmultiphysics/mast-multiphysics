#!/usr/bin/env bash

# Prepare/organize Doxygen generated HTML documentation.
#
#    This script organizes the Doxygen generated HTML documentation and readies it for deployment to
#    a web host.

#    This script should be called during the Travis-CI 'script' phase after documentation is built so that
#    if an error occurs we can catch it on a non-master build (the 'deploy' stage only runs for commits on the
#    `master` branch). We currently call it inside of `ci/build_mast.sh` on the worker that satisfies
#    the condition ${CI_BUILD_DOCS}=true. To actually deploy the documentation (after running this script)
#    call `ci/deploy_docs.sh`.
#

# Organize website content. First make sure we delete any existing website contents (these may exist from
# a previous Travis-CI cache) and then copy generated HTML and other desired files.
rm -rf ${TRAVIS_BUILD_DIR}/website # Don't `|| exit` here. Its not a failure if this directory doesn't exist when
                                   # we try to delete it.
mkdir ${TRAVIS_BUILD_DIR}/website || exit
cd ${TRAVIS_BUILD_DIR}/website || exit
rsync -r ${TRAVIS_BUILD_DIR}/build_rel/doc/doxygen/html/ . || exit
cp ${TRAVIS_BUILD_DIR}/doc/.htaccess . || exit
