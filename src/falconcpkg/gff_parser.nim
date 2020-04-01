# vim: sw=4 ts=4 sts=4 tw=0 et:
import strutils
import streams
import re
from util import nil

let re_attribute = re"cov=(\d+),(\d+),(\d+);cov2=([\d\.]+),([\d\.]+);gaps=(\d+),(\d+)"

type
    GFFCoverageLine* = object
        seq_name*: string
        seq_start*, seq_end*: int
        direction*: string
        cov_min*, cov_max*, n_gaps*, tot_gaps*: int
        cov_median*, cov_mean*, cov_sd*: float

proc parse_gff_line*(line: string): GFFCoverageLine =
    let ld = line.split("\t")
    if ld.len() < 9:
        let msg = "Insufficient number of fields in the GFF line: '{buff}'"
        util.raiseEx(msg)

    var record: GFFCoverageLine

    # Parse the location information.
    record.seq_name = ld[0]
    record.seq_start = parseInt(ld[3]) - 1 # GFF is 1-based.
    record.seq_end = parseInt(ld[4]) # GFF end is inclusive.
    record.direction = ld[6]

    # Parse the coverage.
    let attributes = ld[8]
    var matches: array[7, string]
    let rv_match = match(attributes, re_attribute, matches, 0)
    if not rv_match:
        let msg = "Can't parse coverage from the GFF attributes. Broken file, or file format changed?"
        util.raiseEx(msg)
    record.cov_min = parseInt(matches[0])
    record.cov_median = parseFloat(matches[1])
    record.cov_max = parseInt(matches[2])
    record.cov_mean = parseFloat(matches[3])
    record.cov_sd = parseFloat(matches[4])
    record.n_gaps = parseInt(matches[5])
    record.tot_gaps = parseInt(matches[6])

    return record

iterator yield_gff_record*(sin: Stream): GFFCoverageLine =
    var buff: string
    while(readLine(sin, buff)):
        if buff.len() == 0:
            continue
        if buff[0] == '#':
            continue
        yield(parse_gff_line(buff))

