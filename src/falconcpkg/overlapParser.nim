# vim: sw=4 ts=4 sts=4 tw=0 et:
from os import nil
from strutils import splitWhitespace, parseInt, parseFloat, parseBool
from strformat import fmt
from ./util import isEmptyFile, log
from streams import nil

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

    M4IndexRecord* {.packed.} = object # does not need to be packed unless mem-mapped
                                       #name: string # should be int, and not needed here anyway
        count*: int32
        pos*: int64
        len*: int32
    M4Index* = seq[M4IndexRecord]

proc dumpIndexQuick*(m4idx: M4Index, idx_s: streams.Stream) =
    # Dummy field for Aread.
    for rec in m4idx:
        let desc = "0000 {rec.count} {rec.pos} {rec.len}\n".fmt
        streams.write(idx_s, desc)

proc sscanf(s: cstring, frmt: cstring): cint {.varargs, importc,
        header: "<stdio.h>".}

proc parseM4IndexQuick*(sin: streams.Stream): M4Index =
    # Ignore the Aread field.
    let s_frmt = "%*s %lld %ld %lld"
    var
        line: string
        count: clonglong
        pos: clong
        len: clonglong
        rec: M4IndexRecord
    while streams.readLine(sin, line):
        let scanned = sscanf(line.cstring, s_frmt.cstring,
            addr rec.count, addr rec.pos, addr rec.len)
        if 3 != scanned:
            let msg = "Too few fields ({scanned}) for '{line}'".fmt
            raise newException(util.TooFewFieldsError, msg)
        result.add(rec)

iterator determineNextPile(f: streams.Stream): tuple[count, len: int32] =
    # A "pile" is a block of newline-delimited lines which
    # each start with the same string.
    # On each iteration, set "lenPile" to the length of that string and yield
    # the number of lines in the pile.
    var
        curr: string = "" # Everything starts with this.
        line: string
        index = 0
        count = 0'i32
        lenPile = 0'i32

    while streams.readLine(f, line):
        if not strutils.startsWith(line, curr):
            yield (count, lenPile)
            count = 0
            lenPile = 0
        if 0 == count:
            assert 0 == lenPile
            let iSpace = strutils.find(line, ' ')
            assert iSpace != -1, line
            curr = line[0 .. iSpace]
        lenPile += int32(len(line) + 1) # Note: bug if "\r\l"
        count += 1
    if 0 != count:
        assert 0 != lenPile
        yield (count, lenPile)

proc index*(ovls_s: streams.Stream): M4Index =
    var pos: int64 = 0
    for count, lenPile in determineNextPile(ovls_s):
        let rec = M4IndexRecord(count: count, pos: pos, len: lenPile)
        result.add(rec)
        pos += lenPile

proc getM4Index*(ovls_fn: string): M4Index =
    let idx_fn = ovls_fn & ".idx"
    if os.fileExists(idx_fn):
        #log("Not yet ready to read M4 index '{idx_fn}'. Skipping.".fmt)
        log("Using existing index '{idx_fn}'.".fmt)
        return parseM4IndexQuick(streams.openFilestream(idx_fn, fmRead))
    var ovls_s = streams.openFilestream(ovls_fn)
    let sout = streams.openFilestream(idx_fn, fmWrite)
    if sout.isNil:
        log("Cannot open '{idx_fn}' for write. Not writing index.".fmt)
    result = index(ovls_s)
    if not sout.isNil:
        log("Writing 'quick' index '{idx_fn}', w/ phony Aread names.".fmt)
        dumpIndexQuick(result, sout)

proc parseOverlap*(s: string): Overlap =
    var ovl: Overlap
    let ld = s.splitWhitespace(maxsplit = 12)
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
        ovl.tag = ovl.tagplus.splitWhitespace(maxsplit = 1)[0]
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

iterator getNextPile*(sin: streams.Stream, index: M4Index): seq[Overlap] =
    var
        readA: string
        ov: Overlap
        buff = ""
        ovls = newSeq[Overlap]()

    var totalCount = 0
    #for rec in index:
    #    #echo "count=", rec.count
    #    #for i in 0 ..< count:
    #    #    echo i
    #    totalCount += rec.count
    #echo "totalCount=", totalCount
    #if 0 == len(index):
    #    return
    #echo " index[0]:", index[0]
    assert not sin.isNil
    streams.setPosition(sin, index[0].pos.int)
    for rec in index:
        #echo "  rec:", rec
        assert streams.getPosition(sin) == rec.pos
        #let pile = streams.readStr(sin, rec.len)
        #var ssin = streams.newStringStream(pile)
        for i in 0 ..< rec.count:
            let ok = streams.readLine(sin, buff)
            assert ok, "Failed to read at M4Index {index[i]}".fmt
            ov = parseOverlap(buff)
            ovls.add(ov)
        #if true:
        #    yield @[]
        #    return
        yield ovls
        ovls.setLen(0)

iterator getNextPile*(sin: streams.Stream): seq[Overlap] =
    var
        readA: string
        ov: Overlap
        buff = ""
        ovls = newSeq[Overlap]()

    if streams.readLine(sin, buff):
        try:
            ov = parseOverlap(buff)
            ovls.add(ov)
            readA = ov.Aname
            while(streams.readLine(sin, buff)):
                ov = parseOverlap(buff)
                if ov.Aname != readA:
                    yield ovls
                    readA = ov.Aname
                    ovls.setLen(0)
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
