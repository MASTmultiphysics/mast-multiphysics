#!/usr/bin/env bash

# Deploy organized HTML documentation to the web host.
#
#    This script deploys organized HTML documentation produced by `prepare_docs.sh` to the web host.
#    It is meant to be called during the Travis-CI `deploy` stage and requires a couple of secret
#    tokens to be set in the Travis-CI job. These are the values of WEBSITE_HOST_USER
#    and WEBSITE_HOST_ADDRESS and are specified in the Travis-CI web interface as a secret variable
#    so they are not public in the repo or build logs.

cd ${TRAVIS_BUILD_DIR}/website || exit

# Transfer to main hosting if on master branch and not from a Pull Request into master.
if [ "${TRAVIS_BRANCH}" = "master" ] && [ "${TRAVIS_PULL_REQUEST}" = "false" ]; then
  echo "DEPLOY: Deploying master documentation to https://www.mast-multiphysics.com"
  rsync -e "ssh -o StrictHostKeyChecking=no" -r --delete-after --quiet . ${WEBSITE_HOST_USER}@${WEBSITE_HOST_ADDRESS}:mast-multiphysics.com || exit
fi

# Transfer temporary version of docs for all builds originating within the primary
# MASTmultiphysics/mast-multiphysics repository for bush builds and pull request builds.
# (no forks due to deploy/repo condition in .travis.yml. Storage is identified according to
# the Travis-CI build number and archives are stored for 24 hours before being removed automatically
# from the web server by an external process.
CURRENT_TIME="$(date +%Y_%m_%d-%H:%M:%S:%N)"
echo "${TRAVIS_BUILD_NUMBER}_${CURRENT_TIME}"
