#!/usr/bin/env bash

# Deploy Doxygen generated HTML documentation to Github for hosting.
#
#    This script deploys Doxygen generated HTML documentation to the MAST Github repository for hosting. It is meant
#    to be called during the Travis-CI 'deploy' stage and requires the GH_TOKEN environment variable be set in
#    the Travis-CI job from an account that can push to https://github.com/MASTmultiphysics/MASTmultiphysics.github.io
#
#    Contents of this script were inspired in part by https://gist.github.com/willprice/e07efd73fb7f13f917ea.

# Setup git user for website commit information.
git config --global user.email "travis@travis-ci.org"
git config --global user.name "Travis CI"

# Organize website content. First delete all current contents (except .git directory) and then copy generated
# HTML and other desired files.
mkdir ${TRAVIS_BUILD_DIR}/website
cd ${TRAVIS_BUILD_DIR}/website
rsync -r ${TRAVIS_BUILD_DIR}/build/doc/doxygen/html/ .
cp ${TRAVIS_BUILD_DIR}/doc/README.md .

# Initialize empty git repository. Add all files. Commit.
git init
git add .
git commit --quiet --message "Travis build ${TRAVIS_BUILD_NUMBER}, mast-multiphysics commit ${TRAVIS_COMMIT}"

# Set git remote URL. Force push to GitHub website repository.
# Update git remote URL to include ${GH_TOKEN} key. Push files back to GitHub website repository.
git remote add origin https://${GH_TOKEN}@github.com/MASTmultiphysics/MASTmultiphysics.github.io.git
git push --force origin master

