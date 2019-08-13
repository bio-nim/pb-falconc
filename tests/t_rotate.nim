# vim: sw=4 ts=4 sts=4 tw=0 et:

import falconcpkg/rotate
import unittest
import sets

suite "rotate":
    test "loadWhiteList":
        let got = loadWhiteList("")
        check got == initHashSet[string]()
