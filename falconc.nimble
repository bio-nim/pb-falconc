# Package

version       = "1.4.0"
author        = "Zev Kronenberg"
author        = "Christopher Dunn"
description   = "C Utilities for bioinformatics"
license       = "BSD-3-Clause"
srcDir        = "src"
installDirs   = @["falconcpkg"]
bin           = @["falconc"]


# Dependencies

requires "nim >= 1.2.0", "cligen", "hts", "networkx", "bitvector >= 0.4.10", "msgpack4nim", "threadpools"

task integ, "Runs integration tests":
  var cmd = ""
  cmd = "nim c -r integ-tests/kmer_test1"
  echo cmd
  exec cmd
  #cmd = "nim c -r integ-tests/phasr_test.nim"
  #echo cmd
  #exec cmd
