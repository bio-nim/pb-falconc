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
        bam_filter_clipped(output_fn, input_fn, max_clipping, end_margin, verbose, tags_enrich = true)
        check bam_count(output_fn) == 4
        let osam = readFile(output_fn)
        check osam == enriched_sam
        os.removeFile(output_fn)

    test "bam_filter_clipped_5,-1":
        let
            output_fn = "foo.sam"
            max_clipping = 5
            end_margin = -1
        discard os.tryRemoveFile(output_fn)
        bam_filter_clipped(output_fn, input_fn, max_clipping, end_margin, verbose)
        check bam_count(output_fn) == 4 - 2
        os.removeFile(output_fn)

    test "bam_filter_clipped_6,-1":
        let
            output_fn = "foo.sam"
            max_clipping = 6
            end_margin = -1
        discard os.tryRemoveFile(output_fn)
        bam_filter_clipped(output_fn, input_fn, max_clipping, end_margin, verbose)
        check bam_count(output_fn) == 4
        os.removeFile(output_fn)

    os.removeFile(input_fn)
