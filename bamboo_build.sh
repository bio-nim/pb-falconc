#!/bin/bash
source /mnt/software/Modules/current/init/bash
module load nim/0.20.0
module load gcc ccache

set -vex
make build
