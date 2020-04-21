# vim: sw=4 ts=4 sts=4 tw=0 et:

import falconcpkg/util
import unittest
from sequtils import nil

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
    test "splitWeighted":
        check splitWeighted(0, @[]) == []
        check splitWeighted(0, @[42]) == []
        check splitWeighted(1, @[42]) == [1]
        check splitWeighted(2, @[42, 2]) == [1, 1]
        check splitWeighted(3, @[1, 1, 1]) == [1, 1, 1]
        check splitWeighted(3, @[1, 1, 1, 1]) == [2, 1, 1]
        check splitWeighted(3, @[1, 1]) == [1, 1]
        check splitWeighted(1, @[1, 2, 3, 4]) == [4]
        check splitWeighted(2, @[1, 2, 3, 4]) == [3, 1]
        check splitWeighted(3, @[1, 2, 3, 4]) == [3, 1] # greedy
        check splitWeighted(4, @[1, 2, 3, 4]) == [2, 1, 1] # greedy, so order matters
        check splitWeighted(4, @[4, 3, 2, 1]) == [1, 1, 1, 1] # see?
        check splitWeighted(3, @[4, 3, 2, 1]) == [1, 1, 2]
        check splitWeighted(2, @[4, 3, 2, 1]) == [2, 2]
        check splitWeighted(1, @[4, 3, 2, 1]) == [4]
    test "combineToTarget":
        proc icombineToTarget(t: int, weights: seq[int]): seq[seq[int]] =
            return combineToTarget(t, sequtils.mapIt(weights, int64(it)))
        check icombineToTarget(3, @[2, 2, 2, 2]) == @[@[0, 1], @[2, 3]]
        check icombineToTarget(3, @[2, 2, 2]) == @[@[0, 1], @[2]]
        check icombineToTarget(3, @[2, 2]) == @[@[0, 1]]
        check icombineToTarget(2, @[2, 2, 2, 2]) == @[@[0], @[1], @[2], @[3]]
        check icombineToTarget(1, @[2, 2, 2, 2]) == @[@[0], @[1], @[2], @[3]]
        check icombineToTarget(4, @[2, 2, 2, 2]) == @[@[0, 1], @[2, 3]]
        check icombineToTarget(3, @[1, 2, 3, 4]) == @[@[0, 1], @[2], @[3]]
        check icombineToTarget(3, @[1, 1, 2, 1]) == @[@[0, 1, 2], @[3]]
        check icombineToTarget(3, @[1, 2, 1, 1]) == @[@[0, 1], @[2, 3]]
