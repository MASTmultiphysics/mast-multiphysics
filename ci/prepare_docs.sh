#!/usr/bin/env bash

# Prepare/organize Doxygen generated HTML documentation to be hosted on Github.
#
#    This script organizes the Doxygen generated HTML documentation and commits it to a local .git repository.
#    This script should be called during the Travis-CI 'script' phase after documentation is built so that
#    if an error occurs we can catch it on a non-master build (the 'deploy' stage only runs for commits on the
#    `master` branch). We currently call it inside of `ci/build_mast.sh` on the worker that satisfies
#    the condition ${CI_BUILD_DOCS}=true. To actually deploy the documentation (after running this script)
#    call `ci/deploy_docs.sh`.
#
#    Contents of this script (and `deploy_docs.sh`) were inspired in part by
#    https://gist.github.com/willprice/e07efd73fb7f13f917ea.

# Setup git user for website commit information.
git config --global user.email "travis@travis-ci.org" || exit
git config --global user.name "Travis CI" || exit

# Organize website content. First make sure we delete any existing website contents (these may exist from
# a previous Travis-CI cache) and then copy generated HTML and other desired files.
rm -rf ${TRAVIS_BUILD_DIR}/website # Don't `|| exit` here. Its not a failure if this directory doesn't exist when
                                   # we try to delete it.
mkdir ${TRAVIS_BUILD_DIR}/website || exit
cd ${TRAVIS_BUILD_DIR}/website || exit
rsync -r ${TRAVIS_BUILD_DIR}/build_rel/doc/doxygen/html/ . || exit

# Initialize empty git repository. Add all files. Commit.
git init || exit
git add . || exit
git commit --quiet --message "Travis build ${TRAVIS_BUILD_NUMBER}, mast-multiphysics commit ${TRAVIS_COMMIT}" || exit

# Add remote where we will push website to on actual deployment. GH_TOKEN environment variable is set automatically
# in the Travis-CI, but must from an account that can push to https://github.com/MASTmultiphysics/MASTmultiphysics.github.io
git remote add origin https://${GH_TOKEN}@github.com/MASTmultiphysics/MASTmultiphysics.github.io.git || exit