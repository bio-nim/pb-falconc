from strutils import splitWhitespace, split, parseInt, parseFloat, parseBool
from strformat import fmt
from ./util import isEmptyFile, log
import streams

type
    Overlap* = object
        Aname*: string
        Bname*: string
        score*: int
        idt*: float64
        Arev*: bool
        Astart*: int
        Aend*: int
        Alen*: int
        Brev*: bool
        Bstart*: int
        Bend*: int
        Blen*: int
        tag*: string
        tagplus*: string # tag (redundantly) plus extra columns

proc parseOverlap*(s: string): Overlap =
    var ovl: Overlap
    let ld = s.splitWhitespace(maxsplit=12)
    if ld.len() < 12:
        let msg = "Error parsing ovl (split={ld.len()}): '{s}'".fmt
        log(msg)
    doAssert (ld.len() == 12 or ld.len() == 13)
    ovl.Aname = ld[0]
    ovl.Bname = ld[1]
    ovl.score = parseInt(ld[2])
    ovl.idt = parseFloat(ld[3])
    ovl.Arev = parseBool(ld[4])
    ovl.Astart = parseInt(ld[5])
    ovl.Aend = parseInt(ld[6])
    ovl.Alen = parseInt(ld[7])
    ovl.Brev = parseBool(ld[8])
    ovl.Bstart = parseInt(ld[9])
    ovl.Bend = parseInt(ld[10])
    ovl.Blen = parseInt(ld[11])
    ovl.tag = ""
    ovl.tagplus = ""
    if ld.len() >= 13:
        ovl.tagplus = ld[12]
        ovl.tag = ovl.tagplus.splitWhitespace(maxsplit=1)[0]
    return ovl

proc toString*(o: Overlap): string =
    var strandA = 0
    var strandB = 0
    var tagplus = ""
    if o.Arev: strandA = 1
    if o.Brev: strandB = 1
    if o.tagplus != "": tagplus = " " & o.tagplus
    result = "{o.Aname} {o.Bname} {o.score} {o.idt:.3f} {strandA} {o.Astart} {o.Aend} {o.Alen} {strandB} {o.Bstart} {o.Bend} {o.Blen}{tagplus}".fmt

proc `$`(o: Overlap): string =
    return toString(o)

proc toTuple*(o: Overlap): (string, string, int, float64, bool, int, int, int, bool, int, int, int, string, string) =
    return (o.Aname, o.Bname, o.score, o.idt, o.Arev, o.Astart, o.Aend, o.Alen, o.Brev, o.Bstart, o.Bend, o.Blen, o.tag, o.tagplus)

iterator getNextPile*(sin: Stream): seq[Overlap] =
    var
        readA: string
        ov: Overlap
        buff = ""
        ovls = newSeq[Overlap]()

    if readLine(sin, buff):
        try:
            ov = parseOverlap(buff)
            ovls.add(ov)
            readA = ov.Aname
            while(readLine(sin, buff)):
                ov = parseOverlap(buff)
                if ov.Aname != readA:
                    yield ovls
                    readA = ov.Aname
                    ovls = @[]
                ovls.add(ov)
            yield ovls
        except AssertionError:
            log("Error parsing line: buff='{buff}'".fmt)
            raise

proc parseM4*(sin: streams.Stream): seq[Overlap] =
    for line in streams.lines(sin):
        if strutils.startswith(line, '-'):
            # This can be used to signify EOF, but only if
            # we ever have more trouble with filesystem errors.
            break
        let overlap = parseOverlap(line)
        result.add(overlap)
