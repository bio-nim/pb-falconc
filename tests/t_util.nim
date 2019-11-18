# vim: sw=4 ts=4 sts=4 tw=0 et:

import falconcpkg/util
import unittest

suite "util":
    test "thousands":
        check thousands(0) == "0"
        check thousands(1) == "1"
        check thousands(10) == "10"
        check thousands(100) == "100"
        check thousands(1_000) == "1,000"
        check thousands(10_000) == "10,000"
        check thousands(100_000) == "100,000"
        check thousands(1_000_000) == "1,000,000"
        check thousands(-1_000_000) == "-1,000,000"
        check thousands(-10_000) == "-10,000"
        check thousands(-1_000) == "-1,000"
        check thousands(-1) == "-1"
        check thousands(-0) == "0"
