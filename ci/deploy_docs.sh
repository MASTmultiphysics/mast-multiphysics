#!/usr/bin/env bash

# Deploy organized HTML documentation to the web host.
#
#    This script deploys organized HTML documentation produced by `prepare_docs.sh` to the web host.
#    It is meant to be called during the Travis-CI `deploy` stage and requires a couple of secret
#    tokens to be set in the Travis-CI job. These are the values of WEBSITE_HOST_USER
#    and WEBSITE_HOST_ADDRESS and are specified in the Travis-CI web interface as a secret variable
#    so they are not public in the repo or build logs.

cd ${TRAVIS_BUILD_DIR}/website || exit

# Transfer to main hosting if on master branch.
if [ "${TRAVIS_BRANCH}" = "master" ]; then
  echo "DEPLOY: Deploying master documentation to https://www.mast-multiphysics.com" || exit
  rsync -e "ssh -o StrictHostKeyChecking=no" -r --delete-after --quiet . ${WEBSITE_HOST_USER}@${WEBSITE_HOST_ADDRESS}:mast-multiphysics.com || exit
fi

# Transfer temporary version of docs for all builds originating within the primary
# MASTmultiphysics/mast-multiphysics repository for push builds.
# (no forks due to deploy/repo condition in .travis.yml, and no "pull request" merged builds since
# Travis-CI deploy stage does not run for pull requests to prevent archival of arbitrary artifacts)
# Storage is identified according to the Travis-CI build number and archives are stored for 24 hours before being removed automatically
# from the web server by an external process.
cd .. || exit
CURRENT_TIME="$(date +%Y_%m_%d-%H_%M_%S_%N)"
ARCHIVE_NAME="BUILD_NUM_${TRAVIS_BUILD_NUMBER}_TIMESTAMP_${CURRENT_TIME}.tar.gz"
echo "DEPLOY: Creating website archive: ${ARCHIVE_NAME}" || exit
tar -zcf "${ARCHIVE_NAME}" website || exit
echo "DEPLOY: Transferring archive to https://temp-docs.mast-multiphysics.com" || exit
scp -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null -q -r "${ARCHIVE_NAME}" ${WEBSITE_HOST_USER}@${WEBSITE_HOST_ADDRESS}:temp-docs.mast-multiphysics.com || exit
