# vim: sw=4 ts=4 sts=4 tw=0 et:
import unittest
import "../src/falconcpkg/align"
from os import nil

let sam = """
@SQ	SN:myseq	LN:7
alnA	2064	myseq	1	60	7=	*	0	0	GATTACA	*
alnA	2064	myseq	1	60	5S7=5S	*	0	0	AAAAAGATTACAGGGGG	*
alnA	2064	myseq	1	60	6S7=5S	*	0	0	AAAAAAGATTACAGGGGG	*
alnA	2064	myseq	1	60	5S7=6S	*	0	0	AAAAAGATTACAGGGGGG	*
"""

let enriched_sam = """
@SQ	SN:myseq	LN:7
alnA	2064	myseq	1	60	7=	*	0	0	GATTACA	*	XB:i:0	XE:i:7	XR:i:7	XQ:i:7
alnA	2064	myseq	1	60	5S7=5S	*	0	0	AAAAAGATTACAGGGGG	*	XB:i:0	XE:i:7	XR:i:7	XQ:i:17
alnA	2064	myseq	1	60	6S7=5S	*	0	0	AAAAAAGATTACAGGGGG	*	XB:i:0	XE:i:7	XR:i:7	XQ:i:18
alnA	2064	myseq	1	60	5S7=6S	*	0	0	AAAAAGATTACAGGGGGG	*	XB:i:0	XE:i:7	XR:i:7	XQ:i:18
"""

suite "align":
    let
        input_fn = "input.sam"
        ocount_fn = "outCount.txt"
        n = 4
        verbose = false # true for debugging
    discard os.tryRemoveFile(input_fn)
    writeFile(input_fn, sam)

    test "bam_filter_clipped_5,0":
        let
            output_fn = "foo.sam"
            max_clipping = 5
            end_margin = 0
        discard os.tryRemoveFile(output_fn)
        bam_filter_clipped(ocount_fn, output_fn, input_fn, max_clipping, end_margin,
                Flags_exclude = "0", verbose = verbose, tags_enrich = true)
        check bam_count(output_fn) == n
        let osam = readFile(output_fn)
        check osam == enriched_sam
        os.removeFile(output_fn)

    test "bam_filter_clipped_5,-1":
        let
            output_fn = "foo.sam"
            max_clipping = 5
            end_margin = -1
        discard os.tryRemoveFile(output_fn)
        bam_filter_clipped(ocount_fn, output_fn, input_fn, max_clipping, end_margin,
                Flags_exclude = "0", verbose = verbose)
        check bam_count(output_fn) == n - 2
        os.removeFile(output_fn)

    test "bam_filter_clipped_6,-1":
        let
            output_fn = "foo.sam"
            max_clipping = 6
            end_margin = -1
        discard os.tryRemoveFile(output_fn)
        bam_filter_clipped(ocount_fn, output_fn, input_fn, max_clipping, end_margin,
                Flags_exclude = "0", verbose = verbose)
        check bam_count(output_fn) == n
        os.removeFile(output_fn)

    test "bam_filter_clipped_flags":
        let
            output_fn = "foo.sam"
            max_clipping = 0
            end_margin = 0
        discard os.tryRemoveFile(output_fn)
        bam_filter_clipped(ocount_fn, output_fn, input_fn, max_clipping, end_margin,
                Flags_exclude = "0x800", verbose = verbose, tags_enrich = true)
        check bam_count(output_fn) == 0
        os.removeFile(output_fn)

    test "bam2paf":
        let
            output_p_fn = "foo-p.paf"
            output_a_fn = "foo-a.paf"
        discard os.tryRemoveFile(output_p_fn)
        discard os.tryRemoveFile(output_a_fn)
        bam2paf(input_fn, output_p_fn, output_a_fn)
        let expected_paf = """
alnA	7	0	7	-	myseq	7	0	7	7	7	60
alnA	17	5	12	-	myseq	7	0	7	7	7	60
alnA	18	5	13	-	myseq	7	0	7	7	7	60
alnA	18	6	12	-	myseq	7	0	7	7	7	60
"""
        check open(output_p_fn).readAll() == expected_paf
    os.removeFile(input_fn)

suite "align-utils":
    test "toUint16":
        check toUint16("0xFF") == 255
        check toUint16("0xff") == 255
        check toUint16("0xFFFF") == 65535
        #check toUint16("0x1FFFF") == 65535 # RangeError
        check toUint16("123") == 123
