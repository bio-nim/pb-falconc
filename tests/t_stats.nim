# vim: sw=4 ts=4 sts=4 tw=0 et:
import falconcpkg/stats
import unittest


let expected = """{
    "sum": 120,
    "mean": 40.0,
    "median": 40,
    "max": 50,
    "N100": 30,
    "L100": 3,
    "N90": 30,
    "L90": 3,
    "N50": 40,
    "L50": 2,
    "esize": 42.0
}"""

suite "stats":
    test "calc_stats":
        let reads = @[50'i32, 40, 30]
        let st = stats.calc_stats(reads)

        check stats.to_json(st) == expected
