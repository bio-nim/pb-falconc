# vim: sw=4 ts=4 sts=4 tw=0 et:

import falconcpkg/pbcromwell
import unittest

suite "pbcromwell":
    test "simple":
        check not pbcromwell.should_remove("foo.txt")
        check not pbcromwell.should_remove("foo.las")
        check pbcromwell.should_remove("raw_reads.9.raw_reads.88.las")
        check pbcromwell.should_remove("preads.9.preads.88.las")
    test "dry_run does not raise":
        pbcromwell.remove("foo.las", dry_run=false)
