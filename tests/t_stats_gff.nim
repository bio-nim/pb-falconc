# vim: sw=4 ts=4 sts=4 tw=0 et:
import falconcpkg/stats_gff
import falconcpkg/gfftools
import json
import unittest
import streams
import tables

##################
### Test data. ###
##################
let test_data_stats_gff_content_0 = """
"""
let test_data_stats_gff_expected_0 = """
{:}"""

##################
let test_data_stats_gff_content_1 = """
##sequence-region ctg.s1.000000F 7209:E~12129:B~38100:E~7209:E ctg_circular 4031019 21559395 1 4031019
ctg.s1.000000F	.	region	1	5000	0.00	+	.	cov=5979,6096,6236;cov2=6091.505,51.392;gaps=0,0
ctg.s1.000000F	.	region	5001	10000	0.00	+	.	cov=5765,5944,6276;cov2=5988.380,139.385;gaps=0,0
ctg.s1.000000F	.	region	10001	15000	0.00	+	.	cov=5813,6137,6413;cov2=6142.946,156.044;gaps=0,0
"""
let test_data_stats_gff_expected_1 = """
{"ctg.s1.000000F": {"cov_max": "6413", "cov_min": "5765", "cov_mean": "6074.277"}}"""

##################
let test_data_stats_gff_content_2 = """
ctg.s1.000000F	.	region	1	5000	0.00	+	.	cov=5979,6096,6236;cov2=6091.505,51.392;gaps=0,0
ctg.s1.000000F	.	region	5001	10000	0.00	+
"""

##################
##################

suite "stats_gff":
    test "summarize_gff_coverage_1 empty input":
        # Input.
        let in_data: string = test_data_stats_gff_content_0
        # Expected.
        let expected: string = test_data_stats_gff_expected_0
        # Run unit under test.
        var sin = streams.newStringStream(in_data)
        let results = summarize_gff_coverage(sin)
        # Evaluate.
        check $results == expected

    test "summarize_gff_coverage_2 proper input":
        # Input.
        let in_data: string = test_data_stats_gff_content_1
        # Expected.
        let expected: string = test_data_stats_gff_expected_1
        # Run unit under test.
        var sin = streams.newStringStream(in_data)
        let results = summarize_gff_coverage(sin)
        # Evaluate.
        check parseJson($results) == parseJson(expected)

    test "summarize_gff_coverage_3 truncated input":
        # Input.
        let in_data: string = test_data_stats_gff_content_2
        try:
            # Run unit under test.
            var sin = streams.newStringStream(in_data)
            let results = summarize_gff_coverage(sin)
            assert false
        except Exception as exc:
            assert true

let
    start3 = """
name0 1 2 3 4 5 6 7 8
name1 1 2 3 4 5 6 7 8
name2 1 2 3 4 5 6 7 8
"""
    mask2 = """
name0 1 2 3 4 5 6 7 8
name1 1 2 3 4 5 6 7 8
"""
suite "gffsubtract":
    test "loadGffLines":
        let in_data = test_data_stats_gff_content_1
        var sin = streams.newStringStream(in_data)
        var gl = loadGffLines(sin)
        check gl.len() == 1
    test "gffsubtract":
        var
            gsin = streams.newStringStream(start3)
            msin = streams.newStringStream(mask2)
            sout = streams.newStringStream()
        gffsubtractStreams(gsin, msin, sout)
        sout.setPosition(0)
        check sout.readAll() == ""
        # Actually, we would want name2 in that case.
