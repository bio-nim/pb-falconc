import unittest
import falconcpkg/overlapFilter
from streams import nil

#[
 A-read contained
 000000001 000001050 -9662 98.67 0 0 9656 9656 1 626 10288 12247 contained
]#

suite "overlapFilter parseOK":
    let record1 = "000000001 000001050 -9662 98.67 0 0 9656 9656 1 626 10288 12247 contained"
    let parsed = parseOvl(record1)
    test "ids":
        check parsed.ridA == "000000001"
        check parsed.ridB == "000001050"
    test "score":
        check parsed.score == -9662
    test "idt":
        check parsed.idt == 98.67
    test "strands":
        check parsed.strand1 == false
        check parsed.strand2 == true
    test "starts":
        check parsed.start1 == 0
        check parsed.start2 == 626
    test "ends":
        check parsed.end1 == 9656
        check parsed.end2 == 10288
    test "lengths":
        check parsed.l1 == 9656
        check parsed.l2 == 12247

suite "overlapFilter containedOK":
    test "contained":
        let record1 = "000000001 000001050 -9662 98.67 0 0 9656 9656 1 626 10288 12247 contained"
        let parsed = parseOvl(record1)
        check contained(parsed) == true
    test "notcontained":
        let record1 = "000000001 000001050 -9662 98.67 0 10 9656 9656 1 626 10288 12247 contained"
        let parsed = parseOvl(record1)
        check contained(parsed) == false

suite "overlapFilter idtOK":
    let record1 = "000000001 000001050 -9662 98.67 0 0 9656 9656 1 626 10288 12247 contained"
    let parsed = parseOvl(record1)
    test "ranges":
        check lowIdt(parsed, 0) == false
        check lowIdt(parsed, 100) == true
        check lowIdt(parsed, 98.66) == false
        check lowIdt(parsed, 98.67) == false
        check lowIdt(parsed, 98.68) == true

suite "overlapFilter Fraction":
    test "lOverlap false":
        let record1 = "000000010 000000578 -13654 99.69 0 0 13668 15613 1 0 13654 16307 overlap"
        let parsed = parseOvl(record1)
        check checkFractionOverlap(parsed, 10.0, 90.0) == false
    test "lOverlap true":
        let record1 = "000000010 000000578 -13654 99.69 0 0 15600 15613 1 0 13654 16307 overlap"
        let parsed = parseOvl(record1)
        check checkFractionOverlap(parsed, 10.0, 90.0) == true

suite "overlapFilter terminatedAlignments":
    test "readA pass start":
        let record1 = "000000001 000002664 -5006 98.58 0 0 5004 9656 1 0 5006 10117 overlap"
        let parsed = parseOvl(record1)
        check missingTerminus(parsed) == false
    test "readA pass end":
        let record1 = "000000001 000002664 -5006 98.58 0 10 9656 9656 1 0 5006 10117 overlap"
        let parsed = parseOvl(record1)
        check missingTerminus(parsed) == false
    test "readA fail":
        let record1 = "000000001 000002664 -5006 98.58 0 10 5004 9656 1 0 5006 10117 overlap"
        let parsed = parseOvl(record1)
        check missingTerminus(parsed) == true
    test "readB pass start":
        let record1 = "000000001 000002664 -5006 98.58 0 0 5004 9656 1 0 5006 10117 overlap"
        let parsed = parseOvl(record1)
        check missingTerminus(parsed) == false
    test "readB pass end":
        let record1 = "000000001 000002664 -5006 98.58 0 0 5004 9656 1 10 10117 10117 overlap"
        let parsed = parseOvl(record1)
        check missingTerminus(parsed) == false
    test "readB fail":
        let record1 = "000000001 000002664 -5006 98.58 0 0 5004 9656 1 10 5006 10117 overlap"
        let parsed = parseOvl(record1)
        check missingTerminus(parsed) == true

proc stream(fn: string): auto =
    return streams.newStringStream(readFile(fn))

let fn_105513170 = stream("data/chimera/fn_105513170.ovlp")
let fn_118816930 = stream("data/chimera/fn_118816930.ovlp")
let fn_90965490 = stream("data/chimera/fn_90965490.ovlp")
let fp_17238 = stream("data/chimera/fp_17238.ovlp")
let fp_30455 = stream("data/chimera/fp_30455.ovlp")
let fp_32757 = stream("data/chimera/fp_32757.ovlp")
let tp_5600 = stream("data/chimera/tp_5600.ovlp")
let tp_15210 = stream("data/chimera/tp_15210.ovlp")

suite "chimera filter":
    test "5600 tp marked":
        for i in abOvl(tp_5600):
            check gapInCoverage(i, 4, 90.0) == true
    test "15210 tp marked":
        for i in abOvl(tp_15210):
            check gapInCoverage(i, 4, 90.0) == true
    test "32757 fp not marked":
        for i in abOvl(fp_32757):
            check gapInCoverage(i, 4, 90.0) == false
    test "17238 fp not marked":
        for i in abOvl(fp_17238):
            check gapInCoverage(i, 4, 90.0) == false
    test "30455 FP not marked ":
        for i in abOvl(fp_30455):
            check gapInCoverage(i, 4, 90.0) == false
    test "118816930 CCS tp marked":
        for i in abOvl(fn_118816930):
            check gapInCoverage(i, 4, 90.0) == true
    test "90965490 CCS tp marked":
        for i in abOvl(fn_90965490):
            check gapInCoverage(i, 4, 90.0) == true
    test "105513170 CCS FN, this should be marked, kept for future improvements":
        for i in abOvl(fn_105513170):
            check gapInCoverage(i, 4, 90.0) == false
