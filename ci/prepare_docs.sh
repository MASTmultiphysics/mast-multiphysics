#!/usr/bin/env bash

# Prepare/organized Doxygen generated HTML documentation to Github for hosting by `deploy_docs.sh`.
#
#    This script organizes the Doxygen generated HTML documentation and commits it to local .git repository.
#    Doxygen generated HTML documentation to the MAST Github repository for hosting. This script should be called
#    during the Travis-CI 'script' phase after documentation is built so that if an error occurs we can catch it
#    on a non-master build (the 'deploy' stage only runs for commits on the `master` branch). We currently call it
#    inside of `ci/build_mast.sh` on the worker that satisfies ${CI_BUILD_DOCS}=true.
#
#    Contents of this script (and `deploy_docs.sh`) were inspired in part by
#    https://gist.github.com/willprice/e07efd73fb7f13f917ea.

# Setup git user for website commit information.
git config --global user.email "travis@travis-ci.org" || exit
git config --global user.name "Travis CI" || exit

# Organize website content. First delete all current contents (except .git directory) and then copy generated
# HTML and other desired files.
mkdir ${TRAVIS_BUILD_DIR}/website || exit
cd ${TRAVIS_BUILD_DIR}/website || exit
rsync -r ${TRAVIS_BUILD_DIR}/build_rel/doc/doxygen/html/ . || exit

# Initialize empty git repository. Add all files. Commit.
git init || exit
git add . || exit
git commit --quiet --message "Travis build ${TRAVIS_BUILD_NUMBER}, mast-multiphysics commit ${TRAVIS_COMMIT}" || exit