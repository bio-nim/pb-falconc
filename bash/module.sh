#source /mnt/software/Modules/current/init/bash
module load nim/devel  # https://github.com/nim-lang/Nim/issues/14194
module load gcc ccache

module load htslib  # not really needed until runtime
module load git  # for the revision number
