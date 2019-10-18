#!/bin/bash
type module >& /dev/null || . /mnt/software/Modules/current/init/bash
set -e
source bash/module.sh

set -vex

# Set XDG_CACHE_HOME and CCACHE dirs.

if [[ -z ${bamboo_planRepository_branchName+x} ]]; then
  : #pass
elif [[ ! -d /pbi/flash/bamboo/ccachedir ]]; then
  echo "[WARNING] /pbi/flash/bamboo/ccachedir is missing"
elif [[ $bamboo_planRepository_branchName == develop ]]; then
  export CCACHE_DIR=/pbi/flash/bamboo/ccachedir/${bamboo_shortPlanKey}.${bamboo_shortJobKey}.develop
  export CCACHE_TEMPDIR=/scratch/bamboo.ccache_tempdir
  export XDG_CACHE_HOME=${CCACHE_DIR}
elif [[ $bamboo_planRepository_branchName == master ]]; then
  export CCACHE_DIR=/pbi/flash/bamboo/ccachedir/${bamboo_shortPlanKey}.${bamboo_shortJobKey}.master
  export CCACHE_TEMPDIR=/scratch/bamboo.ccache_tempdir
  export XDG_CACHE_HOME=${CCACHE_DIR}
elif [[ $USER == bamboo ]]; then
  _shortPlanKey=$(echo ${bamboo_shortPlanKey}|sed -e 's/[0-9]*$//')
  export CCACHE_DIR=/pbi/flash/bamboo/ccachedir/${bamboo_shortPlanKey}.${bamboo_shortJobKey}
  export XDG_CACHE_HOME=${CCACHE_DIR}
  if [[ -d /pbi/flash/bamboo/ccachedir/${_shortPlanKey}.${bamboo_shortJobKey}.develop ]]; then
    ( cd /pbi/flash/bamboo/ccachedir/
      cp -a ${_shortPlanKey}.${bamboo_shortJobKey}.develop $CCACHE_DIR
    )
  fi
  export CCACHE_TEMPDIR=/scratch/bamboo.ccache_tempdir
fi

# Test and build

make rsync
make test
#make build  # The develop branch will build and install, so we can skip this.

# Install

case "${bamboo_planRepository_branchName}" in
  develop|master)
    export PREFIX_ARG="/mnt/software/f/falconc/${bamboo_planRepository_branchName}"
    PREFIX=${PREFIX_ARG} make install
    #module load htslib
    ${PREFIX_ARG}/bin/falconc --help
    ;;
  *)
    ;;
esac

#clean removes the results of the test file we need
#git clean -xdf .
