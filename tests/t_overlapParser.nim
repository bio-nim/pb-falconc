import unittest
import falconcpkg/overlapParser as uut
import streams
#from strformat import fmt

let testDataGetNextPile1 = """
15210 10043 -4436 99.910 0 8199 12654 16697 0 8427 12886 12886 u
5600 10005 -639 90.833 0 1944 3264 9722 1 4349 5671 8935 u
5600 10219 -6264 99.936 0 3437 9722 9722 1 245 6528 7272 u
5600 10419 -639 90.833 0 1944 3264 9722 0 1877 3199 8492 u
32757 10030 -88 70.513 0 9470 12284 16890 0 6365 9129 11852 u
32757 10150 -88 70.581 0 9470 12284 16890 1 5052 7819 12414 u
32757 10304 -576 89.165 0 14841 16710 16890 0 2004 3850 4474 u
32757 10882 -822 96.935 0 7946 9088 16890 1 5172 6314 12879 u
"""
let testDataGetNextPile1Expected = [
"""
15210 10043 -4436 99.910 0 8199 12654 16697 0 8427 12886 12886 u
""",
"""
5600 10005 -639 90.833 0 1944 3264 9722 1 4349 5671 8935 u
5600 10219 -6264 99.936 0 3437 9722 9722 1 245 6528 7272 u
5600 10419 -639 90.833 0 1944 3264 9722 0 1877 3199 8492 u
""",
"""
32757 10030 -88 70.513 0 9470 12284 16890 0 6365 9129 11852 u
32757 10150 -88 70.581 0 9470 12284 16890 1 5052 7819 12414 u
32757 10304 -576 89.165 0 14841 16710 16890 0 2004 3850 4474 u
32757 10882 -822 96.935 0 7946 9088 16890 1 5172 6314 12879 u
"""
]


proc helperPileToString(pile: seq[uut.Overlap]): string =
    for ov in pile:
        result.add(uut.toString(ov))
        result.add('\n')

suite "overlapFilter parseOK":
    let record1 = "000000001 000001050 -9662 98.67 0 0 9656 9656 1 626 10288 12247 contained"
    let parsed = uut.parseOverlap(record1)
    test "ids":
        check parsed.Aname == "000000001"
        check parsed.Bname == "000001050"
    test "score":
        check parsed.score == -9662
    test "idt":
        check parsed.idt == 98.67
    test "strands":
        check parsed.Arev == false
        check parsed.Brev == true
    test "starts":
        check parsed.Astart == 0
        check parsed.Bstart == 626
    test "ends":
        check parsed.Aend == 9656
        check parsed.Bend == 10288
    test "lengths":
        check parsed.Alen == 9656
        check parsed.Blen == 12247

suite "overlapParser parseOverlap":
    test "empty input":
        let record = ""
        expect Exception:
            let parsed = uut.parseOverlap(record)

    test "line with 12 fields":
        let record = "000000001 000001050 -9662 98.670 0 0 9656 9656 1 626 10288 12247"
        let parsed = uut.parseOverlap(record)
        check uut.toString(parsed) == record

    test "line with 13 fields":
        let record = "000000001 000001050 -9662 98.670 0 0 9656 9656 1 626 10288 12247 contained"
        let parsed = uut.parseOverlap(record)
        check uut.toString(parsed) == record

    test "wrong number of fields, less than 12":
        let record = "000000001 000001050 -9662 98.670 0 0"
        expect Exception:
            let parsed = uut.parseOverlap(record)

    test "line with extra fields, 14 provided here":
        let record1 = "000000001 000001050 -9662 98.67 0 0 9656 9656 1 626 10288 12247 contained foo"
        let parsed = uut.parseOverlap(record1)
        check parsed.tag == "contained"
        check parsed.tagplus == "contained foo"

    test "line with extra fields, 15 provided here":
        let record1 = "000000001 000001050 -9662 98.67 0 0 9656 9656 1 626 10288 12247 contained foo bar"
        let parsed = uut.parseOverlap(record1)
        check parsed.tag == "contained"
        check parsed.tagplus == "contained foo bar"

suite "overlapParser getNextBlock":
    test "empty input":
        let inStr = ""
        let inStream: Stream = newStringStream(inStr)
        var results: seq[string]
        for pile in uut.getNextPile(inStream):
            results.add(helperPileToString(pile))
        check results.len() == 0

    test "normal input":
        let inStr = testDataGetNextPile1
        let inStream: Stream = newStringStream(inStr)
        var results: seq[string]
        for pile in uut.getNextPile(inStream):
            results.add(helperPileToString(pile))
        check results == testDataGetNextPile1Expected

suite "overlapParser parseM4":
    test "empty input":
        let inStr = ""
        let inStream: Stream = newStringStream(inStr)
        let results = uut.parseM4(inStream)
        check helperPileToString(results) == inStr

    test "normal input":
        let inStr = testDataGetNextPile1
        let inStream: Stream = newStringStream(inStr)
        let results = uut.parseM4(inStream)
        check helperPileToString(results) == inStr

proc checkIndex(m4, expected: string) =
    # Compute M4Index. Then dump. Then
    # parse the dump. So we check both reading and writing.
    let sin = streams.newStringStream(m4)
    var sout = streams.newStringStream()
    let m4idx = index(sin)
    check sin.getPosition() == m4.len
    uut.dumpIndexQuick(m4idx, sout)
    streams.setPosition(sout, 0)
    check streams.readAll(sout) == expected
    let m4idxp = uut.parseM4IndexQuick(streams.newStringStream(expected))
    check m4idx == m4idxp

suite "m4 index":
    test "empty":
        checkIndex("", "")
    test "one":
        let
            m4 = "foo bar 1 2 3\n"
            expected = "0000 1 0 14\n"
        checkIndex(m4, expected)
    test "several":
        let
            m4 = """
001 001 -1 100.000 0 0 0 0 0 0 0 0 overlap
001 002 -1 100.000 0 0 0 0 0 0 0 0 overlap
001 003 -1 100.000 0 0 0 0 0 0 0 0 3
005 008 -1 100.000 0 0 0 0 0 0 0 0 overlap
006 008 -1 100.000 0 0 0 0 0 0 0 0 overlap foo
006 009 -1 100.000 0 0 0 0 0 0 0 0 overlap foo
"""
            expected = """
0000 3 0 123
0000 1 123 43
0000 2 166 94
"""
        checkIndex(m4, expected)
