# vim: sw=4 ts=4 sts=4 tw=0 et:
import falconcpkg/stats
import unittest


let expected_json_on_empty = """{
    "sum": 0,
    "mean": 0.0,
    "median": 0,
    "max": 0,
    "N100": 0,
    "L100": 0,
    "N90": 0,
    "L90": 0,
    "N50": 0,
    "L50": 0,
    "esize": 0.0
}"""
let expected_json = """{
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
let expected_table = """Assembly statistics
       120   sum
        40.0 mean
        40   median
        50   max
        30   min (aka N100)
         3   number (aka L100)
        30   N90
         3   L90
        40   N50
         2   L50
        42.0 E-Size (aka expected size of contig for a random base)
"""

suite "stats":
    test "calc_stats":
        let reads = @[50'i32, 40, 30]
        let st = stats.calc_stats(reads)

        check stats.to_json(st) == expected_json
        check stats.to_table(st) == expected_table
    test "calc_stats empty":
        let reads = newSeq[int32](0)
        let st = stats.calc_stats(reads)

        check stats.to_json(st) == expected_json_on_empty
