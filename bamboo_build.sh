#!/usr/bin/bash
type module >& /dev/null || . /mnt/software/Modules/current/init/bash
set -e
source bash/module.sh

set -vex

make build

case "${bamboo_planRepository_branchName}" in
  develop|master)
    export PREFIX_ARG="/mnt/software/f/falconc/${bamboo_planRepository_branchName}"
    PREFIX=${PREFIX_ARG} make install
    ${PREFIX_ARG}/bin/falconc --help
    ;;
  *)
    ;;
esac
