# vim: sw=4 ts=4 sts=4 tw=0 et:

import falconcpkg/raptor_db

from falconcpkg/util import nil
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
