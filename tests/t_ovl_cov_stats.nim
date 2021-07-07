import unittest
import falconcpkg/ovl_cov_stats
import streams
from strformat import fmt

let testDataCollectFromStream1 = """
15210 10043 -1000 99.900 0 8199 12654 16697 0 8427 12886 12886 u
56000 10005 -1000 99.900 0 3000 10000 10000 1 4000 11000 15000
56000 10005 -1000 99.900 0 5000 10000 10000 0 6000 11000 17000
56000 10005 -1000 99.900 0 0 7000 10000 0 4000 11000 20000
32757 10030 -1000 99.900 0 0 8000 8000 1 4000 12000 15000
32757 10150 -1000 99.900 0 7000 8000 8000 0 0 1000 17000
32757 10304 -1000 99.900 0 0 2000 8000 0 18000 20000 20000
"""
let testDataCollectFromStream1Expected: seq[int] = @[0, 0, 1, 2, 2, 2]

suite "ovl_cov_stats median":
    test "empty input":
        let data: seq[int] = @[]
        let results = ovl_cov_stats.median(data)
        check results == 0.0

    test "odd input length":
        let data: seq[int] = @[1, 2, 3, 4, 5]
        let results = ovl_cov_stats.median(data)
        check results == 3.0

    test "even input length":
        let data: seq[int] = @[1, 2, 3, 4, 5, 6]
        let results = ovl_cov_stats.median(data)
        check results == 3.5

suite "ovl_cov_stats calc_stats":
    test "empty input":
        let data: seq[int] = @[]
        let results = ovl_cov_stats.calc_stats(data)
        let expected = (mean: float(0.0), stddev: float(0.0), median: float(0.0), min: 0, max: 0, n: 0)
        # If this is not converted to string, tuple comparison fails for some reason.
        check $results == $expected

    test "normal input":
        let data: seq[int] = @[10, 10, 11, 12, 1, 15, 15, 15, 18, 12, 11, 20]
        let results = ovl_cov_stats.calc_stats(data)
        let expected = (mean: float(12.5), stddev: float(4.609772228646444), median: float(12.0), min: 1, max: 20, n: 12)
        # If this is not converted to string, tuple comparison fails for some reason.
        check $results == $expected

suite "ovl_cov_stats collect_from_stream":
    # This implicitly tests the count_overlaps proc.

    test "empty stream":
        let inStr = ""
        let inStream: Stream = newStringStream(inStr)
        let results = ovl_cov_stats.collect_from_stream(inStream)
        let expected: seq[int] = @[]
        check results == expected

    test "test data 1":
        let inStr = testDataCollectFromStream1
        let inStream: Stream = newStringStream(inStr)
        let results = ovl_cov_stats.collect_from_stream(inStream)
        let expected = testDataCollectFromStream1Expected
        check results == expected

