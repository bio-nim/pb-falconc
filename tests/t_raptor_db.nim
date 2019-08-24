# vim: sw=4 ts=4 sts=4 tw=0 et:

import falconcpkg/raptor_db

from falconcpkg/util import nil
import json
import streams
import strutils
import unittest

let content0 = """\
V	0.2
F	0	test-data/raptordb-fetch/test-2.block.0	fasta
S	0	m1/3005/0_5852	5852	0	0	5929
S	1	m1/3414/0_11983	11983	0	5929	12061
S	2	m1/3820/0_24292	24292	0	17990	24370
S	3	m1/3981/0_5105	5105	0	42360	5182
B	0	0	3	47232
"""

let content1 = """\
V	0.2
F	0	test-data/raptordb-fetch/test-2.block.0	fasta
S	1	$#	11	22	33	44
B	0	0	1	11
"""

let content2 = """\
V	0.2
F	0	test-data/raptordb-fetch/test-2.block.0	fasta
S	1	header	11	22	33
B	0	0	1	11
"""

# Total length sum is 550000 bp.
# The sequences will be sorted by size in descending order internally.
let input_data_str = """\
V	0.2
F	0	subreads.bam	bam
S	0	seq/0/0_10000	10000	0	0	10000
S	1	seq/1/0_20000	20000	0	10000	20000
S	2	seq/2/0_30000	30000	0	20000	30000
S	3	seq/3/0_40000	40000	0	30000	40000
S	4	seq/4/0_50000	50000	0	40000	50000
S	5	seq/5/0_60000	60000	0	50000	60000
S	6	seq/6/0_70000	70000	0	60000	70000
S	7	seq/7/0_80000	80000	0	70000	80000
S	8	seq/8/0_90000	90000	0	80000	90000
S	9	seq/9/0_100000	100000	0	90000	100000
B	0	0	10	100000
"""

let test_data = %* [
    { # 9x coverage of a 10kbp genome is 90kbp, and the longest read is 100kbp.
        "input_data": input_data_str,
        "genome_size": 10000,
        "coverage": 9,
        "cmd_params": ["--fail-low-cov"],
        "exp_success": true,
        "exp_out": "100000",
    },
    { # 15x coverage of a 10kbp genome is 150kbp.
        "input_data": input_data_str,
        "genome_size": 10000,
        "coverage": 15,
        "cmd_params": ["--fail-low-cov"],
        "exp_success": true,
        "exp_out": "90000",
    },
    { # 40x coverage of a 10kbp genome is 400kbp.
        "input_data": input_data_str,
        "genome_size": 10000,
        "coverage": 40,
        "cmd_params": ["--fail-low-cov"],
        "exp_success": true,
        "exp_out": "60000",
    },
    { # 40x coverage of a 10kbp genome is 400kbp. Not using "--fail-low-cov".
        "input_data": input_data_str,
        "genome_size": 10000,
        "coverage": 40,
        "cmd_params": [],
        "exp_success": true,
        "exp_out": "60000",
    },
    { # 60x coverage of a 10kbp genome is 600kbp. Expected failure.
        "input_data": input_data_str,
        "genome_size": 10000,
        "coverage": 60,
        "cmd_params": ["--fail-low-cov"],
        "exp_success": false,
        "exp_out": "10000",
    },
    { # 60x coverage of a 10kbp genome is 600kbp. Should not fail, because the "--fail-low-cov" is not specified.
        "input_data": input_data_str,
        "genome_size": 10000,
        "coverage": 60,
        "cmd_params": [],
        "exp_success": true,
        "exp_out": "10000",
    },
    { # Zero cutoff. Expected failure.
        "input_data": input_data_str,
        "genome_size": 10000,
        "coverage": 0,
        "cmd_params": ["--fail-low-cov"],
        "exp_success": false,
        "exp_out": "0",
    },
    { # Negative cutoff. Expected failure.
        "input_data": input_data_str,
        "genome_size": 10000,
        "coverage": -100,
        "cmd_params": ["--fail-low-cov"],
        "exp_success": false,
        "exp_out": "0",
    },
    { # Zero genome size. Expected failure.
        "input_data": input_data_str,
        "genome_size": 0,
        "coverage": 30,
        "cmd_params": ["--fail-low-cov"],
        "exp_success": false,
        "exp_out": "0",
    },
    { # Negative genome size. Expected failure.
        "input_data": input_data_str,
        "genome_size": -100,
        "coverage": 30,
        "cmd_params": ["--fail-low-cov"],
        "exp_success": false,
        "exp_out": "0",
    },
]

suite "raptor_db":

    test "load_db":
        let sin = streams.newStringStream(content0)
        let reqs = load_rdb(sin)
        let n = len(reqs)
        assert 4 == n
        assert "m1/3005/0_5852" == reqs[0].header
    test "get_length_cutoff":
        let sin = streams.newStringStream(content0)

        block:
            sin.setPosition(0)
            expect(AssertionError):
                discard get_length_cutoff(sin, 0, 30)
        block:
            sin.setPosition(0)
            expect(AssertionError):
                discard get_length_cutoff(sin, 10_000, 0)

        block:
            sin.setPosition(0)
            let cutoff = get_length_cutoff(sin, 36_000, 1)
            assert 11983 == cutoff

        block:
            sin.setPosition(0)
            expect(util.GenomeCoverageError):
                discard get_length_cutoff(sin, 10_000_000, 30,
                        fail_low_cov = true)

        block:
            sin.setPosition(0)
            let cutoff = get_length_cutoff(sin, 10_000_000, 30,
                    fail_low_cov = false)
            assert 5105 == cutoff
    test "get_length_cutoff_long_header":
        block:
            let ok_header = strutils.repeat("A", 1023)
            let content = content1 % [ok_header]
            let sin = streams.newStringStream(content)

            discard get_length_cutoff(sin, 100, 30)
        block:
            let long_header = strutils.repeat("A", 1024)
            let content = content1 % [long_header]
            let sin = streams.newStringStream(content)

            expect(util.FieldTooLongError):
                discard get_length_cutoff(sin, 100, 30)
    test "get_length_cutoff_missing_field":
        let sin = streams.newStringStream(content2)

        expect(util.TooFewFieldsError):
            discard get_length_cutoff(sin, 100, 30)
