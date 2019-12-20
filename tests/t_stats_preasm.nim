# vim: sw=4 ts=4 sts=4 tw=0 et:
import falconcpkg/stats_preasm
from falconcpkg/stats_preasm import SeqLenTuple
import unittest
import streams

#######################################################
### Test data for "input_parse_seq_len_tuples_1".   ###
#######################################################
# Normal test case.
let input_parse_seq_len_tuples_1 = """
V	0.2
F	0	subreads.bam	bam
S	0	seq/0/0_10000	10000	0	0	10000
S	1	seq/1/0_20000	20000	0	10000	20000
S	2	seq/2/0_30000	30000	0	30000	30000
B	0	0	3	60000
"""
let expected_parse_seq_len_tuples_1: seq[SeqLenTuple] = @[
    (0.int64, "seq/0/0_10000", 10000.int32),
    (1.int64, "seq/1/0_20000", 20000.int32),
    (2.int64, "seq/2/0_30000", 30000.int32),
]

# Malformed sequence line, has spaces instead of tabs.
let input_parse_seq_len_tuples_2 = """
S   0   seq/0/0_10000   10000   0   0   10000
"""

# Empty input.
let input_parse_seq_len_tuples_3 = """
"""
let expected_parse_seq_len_tuples_3: seq[SeqLenTuple] = @[]

#########################################
### Test data for "cutoff_lengths".   ###
#########################################
# Empty input.
let input_cutoff_lengths_1_cutoff: int32 = 20000
let input_cutoff_lengths_1_lengths: seq[int32] = @[]
let input_cutoff_lengths_1_expected: seq[int32] = @[]

# Empty input.
let input_cutoff_lengths_2_cutoff: int32 = 20000
let input_cutoff_lengths_2_lengths: seq[int32] = @[1000.int32, 20000.int32, 12345.int32, 1234567.int32, 1000000.int32]
let input_cutoff_lengths_2_expected: seq[int32] = @[20000.int32, 1234567.int32, 1000000.int32]

#########################################
### Test data for "run_from_streams". ###
#########################################
# Normal test case, all seed reads are error corrected.
let input_run_from_streams_1_rawreads = """
V	0.2
F	0	subreads.bam	bam
S	0	seq/0/0_10000	10000	0	0	10000
S	1	seq/1/0_20000	20000	0	10000	20000
S	2	seq/2/0_30000	30000	0	30000	30000
B	0	0	3	60000
"""
let input_run_from_streams_1_preads = """
V	0.2
F	0	subreads.bam	bam
S	0	0	9000	0	0	9000
S	1	1	18000	0	10000	18000
S	2	2	25000	0	27000	25000
B	0	0	3	52000
"""
let input_run_from_streams_1_genome_length = 52000
let input_run_from_streams_1_cutoff_length = 0
let expected_run_from_streams_1 = """{
    "genome_length": 52000,
    "length_cutoff": 0,
    "raw_reads": 3,
    "raw_bases": 60000,
    "raw_mean": 20000.0,
    "raw_n50": 30000,
    "raw_p95": 30000,
    "raw_coverage": 1.15,
    "raw_esize": 23333.33,
    "seed_reads": 3,
    "seed_bases": 60000,
    "seed_mean": 20000.0,
    "seed_n50": 30000,
    "seed_p95": 30000,
    "seed_coverage": 1.15,
    "seed_esize": 23333.33,
    "preassembled_reads": 3,
    "preassembled_bases": 52000,
    "preassembled_mean": 17333.33,
    "preassembled_n50": 18000,
    "preassembled_p95": 25000,
    "preassembled_coverage": 1.0,
    "preassembled_esize": 19807.69,
    "preassembled_yield": 0.87,
    "preassembled_seed_fragmentation": 1.0,
    "preassembled_seed_truncation": 2666.67
}"""

# Normal test case, one seed read was dropped.
let input_run_from_streams_2_rawreads = input_run_from_streams_1_rawreads
let input_run_from_streams_2_preads = """
V	0.2
F	0	subreads.bam	bam
S	1	1	18000	0	10000	18000
S	2	2	25000	0	27000	25000
B	0	0	2	43000
"""
let input_run_from_streams_2_genome_length = 52000
let input_run_from_streams_2_cutoff_length = 0
let expected_run_from_streams_2 = """{
    "genome_length": 52000,
    "length_cutoff": 0,
    "raw_reads": 3,
    "raw_bases": 60000,
    "raw_mean": 20000.0,
    "raw_n50": 30000,
    "raw_p95": 30000,
    "raw_coverage": 1.15,
    "raw_esize": 23333.33,
    "seed_reads": 3,
    "seed_bases": 60000,
    "seed_mean": 20000.0,
    "seed_n50": 30000,
    "seed_p95": 30000,
    "seed_coverage": 1.15,
    "seed_esize": 23333.33,
    "preassembled_reads": 2,
    "preassembled_bases": 43000,
    "preassembled_mean": 21500.0,
    "preassembled_n50": 25000,
    "preassembled_p95": 25000,
    "preassembled_coverage": 0.83,
    "preassembled_esize": 22069.77,
    "preassembled_yield": 0.72,
    "preassembled_seed_fragmentation": 1.0,
    "preassembled_seed_truncation": 3500.0
}"""

# Empty inputs.
let input_run_from_streams_3_rawreads = ""
let input_run_from_streams_3_preads = "" 
let input_run_from_streams_3_genome_length = 52000
let input_run_from_streams_3_cutoff_length = 0
let expected_run_from_streams_3 = """{
    "genome_length": 52000,
    "length_cutoff": 0,
    "raw_reads": 0,
    "raw_bases": 0,
    "raw_mean": 0.0,
    "raw_n50": 0,
    "raw_p95": 0,
    "raw_coverage": 0.0,
    "raw_esize": 0.0,
    "seed_reads": 0,
    "seed_bases": 0,
    "seed_mean": 0.0,
    "seed_n50": 0,
    "seed_p95": 0,
    "seed_coverage": 0.0,
    "seed_esize": 0.0,
    "preassembled_reads": 0,
    "preassembled_bases": 0,
    "preassembled_mean": 0.0,
    "preassembled_n50": 0,
    "preassembled_p95": 0,
    "preassembled_coverage": 0.0,
    "preassembled_esize": 0.0,
    "preassembled_yield": 0.0,
    "preassembled_seed_fragmentation": 1.0,
    "preassembled_seed_truncation": -1.0
}"""

# Cutoff value is larger than any sequence. There should be 0 seed reads and 0 preads.
let input_run_from_streams_4_rawreads = input_run_from_streams_1_rawreads
let input_run_from_streams_4_preads = """
"""
let input_run_from_streams_4_genome_length = 52000
let input_run_from_streams_4_cutoff_length = 90000
let expected_run_from_streams_4 = """{
    "genome_length": 52000,
    "length_cutoff": 90000,
    "raw_reads": 3,
    "raw_bases": 60000,
    "raw_mean": 20000.0,
    "raw_n50": 30000,
    "raw_p95": 30000,
    "raw_coverage": 1.15,
    "raw_esize": 23333.33,
    "seed_reads": 0,
    "seed_bases": 0,
    "seed_mean": 0.0,
    "seed_n50": 0,
    "seed_p95": 0,
    "seed_coverage": 0.0,
    "seed_esize": 0.0,
    "preassembled_reads": 0,
    "preassembled_bases": 0,
    "preassembled_mean": 0.0,
    "preassembled_n50": 0,
    "preassembled_p95": 0,
    "preassembled_coverage": 0.0,
    "preassembled_esize": 0.0,
    "preassembled_yield": 0.0,
    "preassembled_seed_fragmentation": 1.0,
    "preassembled_seed_truncation": -1.0
}"""

suite "stats_preasm__parse_seq_len_tuples":
    test "parse_seq_len_tuples_1":
        # Normal test case.
        let input = input_parse_seq_len_tuples_1
        let expected = expected_parse_seq_len_tuples_1
        var sin = newStringStream(input)
        defer: sin.close()
        # Run the test.
        let result = stats_preasm.parse_seq_len_tuples(sin)
        check result == expected

    test "parse_seq_len_tuples_2":
        # Exception should be thrown, the input is ill formatted.
        let input = input_parse_seq_len_tuples_2
        var sin = newStringStream(input)
        defer: sin.close()
        # Run the test.
        expect Exception:
            let result = stats_preasm.parse_seq_len_tuples(sin)

    test "parse_seq_len_tuples_3":
        # Empty input
        let input = input_parse_seq_len_tuples_3
        let expected = expected_parse_seq_len_tuples_3
        var sin = newStringStream(input)
        defer: sin.close()
        # Run the test.
        let result = stats_preasm.parse_seq_len_tuples(sin)
        check result == expected

suite "stats_preasm__cutoff_lengths":
    test "parse_seq_len_tuples_1":
        # Empty input.
        let input_cutoff = input_cutoff_lengths_1_cutoff
        let input_lengths = input_cutoff_lengths_1_lengths
        let expected = input_cutoff_lengths_1_expected
        # Run the test.
        let result = stats_preasm.cutoff_lengths(input_cutoff, input_lengths)
        check result == expected
    test "parse_seq_len_tuples_2":
        # Normal test case.
        let input_cutoff = input_cutoff_lengths_2_cutoff
        let input_lengths = input_cutoff_lengths_2_lengths
        let expected = input_cutoff_lengths_2_expected
        # Run the test.
        let result = stats_preasm.cutoff_lengths(input_cutoff, input_lengths)
        check result == expected

suite "stats_preasm":
    test "run_from_streams_1":
        # Normal test case. 3 input raw reads, seed length cutoff is 0 so there
        # are in total 3 seed reads, and all 3 were error corrected.
        let input_rawreads_db = input_run_from_streams_1_rawreads
        let input_preads_db = input_run_from_streams_1_preads
        let input_genome_length = input_run_from_streams_1_genome_length.int32
        let input_cutoff_length = input_run_from_streams_1_cutoff_length.int32
        let expected = expected_run_from_streams_1
        # Initialize the stream for rawreads.
        let stream_rawreads: Stream = newStringStream(input_rawreads_db)
        defer: stream_rawreads.close()
        # Initialize the stream for preads.
        let stream_preads: Stream = newStringStream(input_preads_db)
        defer: stream_preads.close()
        # Run the test.
        let result = stats_preasm.run_from_streams(stream_rawreads, stream_preads, input_genome_length, input_cutoff_length)
        check result == expected

    test "run_from_streams_2":
        # Same as the previous test (3 raw reads, 3 seed reads), but only 2 preads.
        let input_rawreads_db = input_run_from_streams_2_rawreads
        let input_preads_db = input_run_from_streams_2_preads
        let input_genome_length = input_run_from_streams_2_genome_length.int32
        let input_cutoff_length = input_run_from_streams_2_cutoff_length.int32
        let expected = expected_run_from_streams_2
        # Initialize the stream for rawreads.
        let stream_rawreads: Stream = newStringStream(input_rawreads_db)
        defer: stream_rawreads.close()
        # Initialize the stream for preads.
        let stream_preads: Stream = newStringStream(input_preads_db)
        defer: stream_preads.close()
        # Run the test.
        let result = stats_preasm.run_from_streams(stream_rawreads, stream_preads, input_genome_length, input_cutoff_length)
        check result == expected

    test "run_from_streams_3":
        # Empty inputs.
        let input_rawreads_db = input_run_from_streams_3_rawreads
        let input_preads_db = input_run_from_streams_3_preads
        let input_genome_length = input_run_from_streams_3_genome_length.int32
        let input_cutoff_length = input_run_from_streams_3_cutoff_length.int32
        let expected = expected_run_from_streams_3
        # Initialize the stream for rawreads.
        let stream_rawreads: Stream = newStringStream(input_rawreads_db)
        defer: stream_rawreads.close()
        # Initialize the stream for preads.
        let stream_preads: Stream = newStringStream(input_preads_db)
        defer: stream_preads.close()
        # Run the test.
        let result = stats_preasm.run_from_streams(stream_rawreads, stream_preads, input_genome_length, input_cutoff_length)
        check result == expected

    test "run_from_streams_4":
        # Large seed length cutoff, larger than any raw read.
        let input_rawreads_db = input_run_from_streams_4_rawreads
        let input_preads_db = input_run_from_streams_4_preads
        let input_genome_length = input_run_from_streams_4_genome_length.int32
        let input_cutoff_length = input_run_from_streams_4_cutoff_length.int32
        let expected = expected_run_from_streams_4
        # Initialize the stream for rawreads.
        let stream_rawreads: Stream = newStringStream(input_rawreads_db)
        defer: stream_rawreads.close()
        # Initialize the stream for preads.
        let stream_preads: Stream = newStringStream(input_preads_db)
        defer: stream_preads.close()
        # Run the test.
        let result = stats_preasm.run_from_streams(stream_rawreads, stream_preads, input_genome_length, input_cutoff_length)
        check result == expected
