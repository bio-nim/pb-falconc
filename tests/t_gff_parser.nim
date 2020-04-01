# vim: sw=4 ts=4 sts=4 tw=0 et:
import falconcpkg/gff_parser
import unittest
import streams

##################
### Test data. ###
##################
let test_data_gff_parser_content_0 = """

##sequence-region ctg.s1.000000F 7209:E~12129:B~38100:E~7209:E ctg_circular 4031019 21559395 1 4031019
# This is a comment!
"""

let test_data_gff_parser_content_1 = """
# This is a comment!
ctg.s1.000000F	.	region	1	5000	0.00	+
"""

let test_data_gff_parser_content_2 = """
##sequence-region ctg.s1.000000F 7209:E~12129:B~38100:E~7209:E ctg_circular 4031019 21559395 1 4031019
ctg.s1.000000F	.	region	1	5000	0.00	+	.	cov=5979,6096,6236;cov2=6091.505,51.392;gaps=0,0
ctg.s1.000000F	.	region	5001	10000	0.00	+	.	cov=5765,5944,6276;cov2=5988.380,139.385;gaps=0,0
ctg.s1.000000F	.	region	10001	15000	0.00	+	.	cov=5813,6137,6413;cov2=6142.946,156.044;gaps=0,0
"""
    ##################
    ##################

suite "gff_parser":
    test "parse_line_1 empty input":
        # Empty input.
        let in_data = ""
        let expected: GFFCoverageLine = GFFCoverageLine()
        try:
            let result = parse_gff_line(in_data)
            assert false
        except Exception as exc:
            assert true

    test "parse_line_2 proper input":
        # Proper input.
        let in_data = "ctg.s1.000000F	.	region	1785001	1790000	0.00	+	.	cov=4451,4564,4718;cov2=4567.164,49.624;gaps=0,0"
        let expected: GFFCoverageLine = GFFCoverageLine(
            seq_name: "ctg.s1.000000F",
            seq_start: 1785000,
            seq_end: 1790000,
            direction: "+",
            cov_min: 4451,
            cov_median: 4564,
            cov_max: 4718,
            cov_mean: 4567.164,
            cov_sd: 49.624,
            n_gaps: 0,
            tot_gaps: 0
        )
        let result = parse_gff_line(in_data)

    test "parse_line_3 wrong number of fields":
        # Empty input.
        let in_data = "ctg.s1.000000F	.	region	1785001	1790000	0.00"
        let expected: GFFCoverageLine = GFFCoverageLine()
        try:
            let result = parse_gff_line(in_data)
            assert false
        except Exception as exc:
            assert true

    test "yield_gff_record_1 empty input":
        # Input.
        let in_data = ""
        # Expected.
        let expected_num_records = 0
        # Run unit under test and collect all records.
        var records = newSeq[GFFCoverageLine]()
        var sin = streams.newStringStream(in_data)
        for record in yield_gff_record(sin):
            records.add(record)
        # Evaluate.
        check records.len() == expected_num_records

    test "yield_gff_record_2 only comment lines, zero records should be read":
        # Input.
        let in_data = test_data_gff_parser_content_0
        # Expected.
        let expected_num_records = 0
        # Run unit under test and collect all records.
        var records = newSeq[GFFCoverageLine]()
        var sin = streams.newStringStream(in_data)
        for record in yield_gff_record(sin):
            records.add(record)
        # Evaluate.
        check records.len() == expected_num_records


    test "yield_gff_record_3 truncated file, this should raise":
        # Input.
        let in_data = test_data_gff_parser_content_1
        # Run unit under test and collect all records.
        var records = newSeq[GFFCoverageLine]()
        var sin = streams.newStringStream(in_data)
        try:
            for record in yield_gff_record(sin):
                records.add(record)
        except Exception as exc:
            assert true

    test "yield_gff_record_2 proper input":
        # Input.
        let in_data = test_data_gff_parser_content_2

        # Expected.
        let expected_records = @[
            GFFCoverageLine(
                seq_name: "ctg.s1.000000F",
                seq_start: 0,
                seq_end: 5000,
                direction: "+",
                cov_min: 5979,
                cov_median: 6096,
                cov_max: 6236,
                cov_mean: 6091.505,
                cov_sd: 51.392,
                n_gaps: 0,
                tot_gaps: 0
            ),
            GFFCoverageLine(
                seq_name: "ctg.s1.000000F",
                seq_start: 5000,
                seq_end: 10000,
                direction: "+",
                cov_min: 5765,
                cov_median: 5944,
                cov_max: 6276,
                cov_mean: 5988.380,
                cov_sd: 139.385,
                n_gaps: 0,
                tot_gaps: 0
            ),
            GFFCoverageLine(
                seq_name: "ctg.s1.000000F",
                seq_start: 10000,
                seq_end: 15000,
                direction: "+",
                cov_min: 5813,
                cov_median: 6137,
                cov_max: 6413,
                cov_mean: 6142.946,
                cov_sd: 156.044,
                n_gaps: 0,
                tot_gaps: 0
            )
        ]

        # Run unit under test and collect all records.
        var records = newSeq[GFFCoverageLine]()
        var sin = streams.newStringStream(in_data)
        for record in yield_gff_record(sin):
            records.add(record)

        # Evaluate.
        check records.len() == expected_records.len()
        for i in 0..<records.len():
            check records[i].seq_name == expected_records[i].seq_name
            check records[i].seq_start == expected_records[i].seq_start
            check records[i].seq_end == expected_records[i].seq_end
            check records[i].direction == expected_records[i].direction
            check records[i].cov_min == expected_records[i].cov_min
            check records[i].cov_median == expected_records[i].cov_median
            check records[i].cov_max == expected_records[i].cov_max
            check records[i].cov_mean == expected_records[i].cov_mean
            check records[i].cov_sd == expected_records[i].cov_sd
            check records[i].n_gaps == expected_records[i].n_gaps
            check records[i].tot_gaps == expected_records[i].tot_gaps

