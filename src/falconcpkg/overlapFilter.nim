import sets
import tables
from strutils import split, parseInt, parseFloat, parseBool
from sequtils import keepIf
from strformat import fmt
from ./util import isEmptyFile
import msgpack4nim
import streams
import json
import osproc
import algorithm


#[
 TODOs:
 1. Allow non-symmetrical overlaps
 2. Make sure the string buffer is an okay length.
 3. Maybe write n-best overlaps code
]#

#[
LA4falcon Columns:
 1. readId1
 2. readId2
 3. score?
 4. % similarity -- estimated
 5. strand 1
 6. startRead1
 7. endRead1
 8. lenght read1
 9. strand 2
 10. startRead2
 11. endRead2
 12. length read2
 13. tag
]#

#Bitfield for read classification
const CREAD = 1 # contained read
const SREAD = 2 # small read
const ZREAD = 4 # self hit
const LREAD = 8 # one side has a low edge count
const HREAD = 16 # one side has a high edge count
const BREAD = 32 # read balance is off
const GREAD = 64 # read has unspanned bases, putative chimera

type
    overlap* = object
        ridA*: string
        ridB*: string
        score*: int
        idt*: float
        strand1*: bool
        start1*: int
        end1*: int
        l1*: int
        strand2*: bool
        start2*: int
        end2*: int
        l2*: int
        tag: string

    coverageInfo = tuple[gap: bool, lowCov: bool, highCov: bool,
            balance: bool]

    ScoredOverlap = tuple
        ovlpLen: int
        mRange: int
        o: overlap



#var readsToFilter = initTable[string, int]()

proc parseOvl*(s: string): overlap =
    var ovl: overlap
    let ld = s.split(" ")
    doAssert ld.len() == 13
    ovl.ridA = ld[0]
    ovl.ridB = ld[1]
    ovl.score = parseInt(ld[2])
    ovl.idt = parseFloat(ld[3])
    ovl.strand1 = parseBool(ld[4])
    ovl.start1 = parseInt(ld[5])
    ovl.end1 = parseInt(ld[6])
    ovl.l1 = parseInt(ld[7])
    ovl.strand2 = parseBool(ld[8])
    ovl.start2 = parseInt(ld[9])
    ovl.end2 = parseInt(ld[10])
    ovl.l2 = parseInt(ld[11])
    ovl.tag = ld[12]
    return ovl

proc `$`(o: overlap): string =
    var strand1 = 0
    var strand2 = 0
    if o.strand1: strand1 = 1
    if o.strand2: strand2 = 1

    result = "{o.ridA} {o.ridB} {o.score} {o.idt:.3f} {strand1} {o.start1} {o.end1} {o.l1} {strand2} {o.start2} {o.end2} {o.l2} {o.tag}".fmt

proc lowIdt*(o: overlap, idt: float): bool =
    ##Check if percent identity is below threshold
    #This is a overlap filter
    if o.idt < idt:
        return true
    return false

proc checkFractionOverlap*(o: overlap, lowOverlap: float,
        highOverlap: float): bool =
    ##Check for percentage of two reads that overlap
    #This is an overlap filter - new as of June 2019
    let perA = (float(o.end1 - o.start1) / (float(o.l1))) * 100
    let perB = (float(o.end2 - o.start2) / (float(o.l2))) * 100
    if (perA < lowOverlap) or (perA > highOverlap) or (perB < lowOverlap) or (
            perB > highOverlap):
        return true
    return false

proc small*(o: overlap, l: int): bool =
    ##Tests if a read is too small
    #This is a read filter
    if o.l1 < l:
        return true
    else:
        return false

proc contained*(o: overlap): bool =
    ##Tests if A-read is contained within B
    #This is a read filter
    if o.start1 > 0:
        return false
    if o.end1 != o.l1:
        return false
    else:
        return true

proc missingTerminus*(o: overlap): bool =
    ##Test that overlap reaches the end of both reads
    #This is an overlap filter
    if o.start1 != 0 and o.end1 != o.l1:
        return true
    if o.start2 != 0 and o.end2 != o.l2:
        return true
    else:
        return false



iterator abOvl(): seq[overlap] =
    var readA: string
    #The buffer len of 500 should be adjusted if more fields are added.
    var buff = newString(500)
    var ovls = newSeq[overlap]()
    discard readLine(stdin, buff)
    var ov = parseOvl(buff)
    readA = ov.ridA
    ovls.add(ov)
    while(readLine(stdin, buff)):
        ov = parseOvl(buff)
        if ov.ridA != readA:
            yield ovls
        if ov.ridA != readA:
            readA = ov.ridA
            ovls = @[]
        ovls.add(ov)
    yield ovls

proc summarize(filterLog: string, readsToFilter: var Table[string, int]) =
    var readFilterCounts = initCountTable[string]()
    for k, v in readsToFilter:
        if (v and CREAD) > 0: readFilterCounts.inc("contained read")
        if (v and SREAD) > 0: readFilterCounts.inc("small read")
        if (v and LREAD) > 0: readFilterCounts.inc("low overlap count")
        if (v and HREAD) > 0: readFilterCounts.inc("high overlap count")
        if (v and BREAD) > 0: readFilterCounts.inc("not balanced")
        if (v and GREAD) > 0: readFilterCounts.inc("gap in coverage")

    var fout = open(filterLog, fmWrite)
    fout.writeLine("filtered reads:", len(readsToFilter))
    for k, v in readFilterCounts:
        fout.writeLine("flag: ", k, " :count:", v)
    var report = """



    Small read:
     The read length is smaller than the default or user defined
     cutoff (--minLen 6,000bp)

    Contained read:
     In the example below A-read is contained within B-read.
     If the B-read has been marked as a repeat or too small the A-read is no
     longer contained within B.

     -----------B-----------
       ---A---

    Low/High overlap count read:
     The number of overlaping B-reads is low/high on either side of the A-read.
     This filter is controlled by --minCov/--maxCov with defaults 2/200. In the
     example below the read fails default filters because the 3prime of A-read
     only has one B-read covering it. The B-read alignments must be dovetail.

     ---B---              ---B---
      -----------A-----------
     ---B---

    Not balanced read:
     The absolute difference between the number of 5prime and 3prime overlaps
     counts are greater than the default or user defined cutoff (--maxDiff 100).

    Gap in coverage read:
     There is a region along the A-read that has low coverage. The --gapFilt
     flag toggles this filter. The cutoff is controlled by (--minDepth 2).
     Chimeric reads tend to have gaps in coverage since they contain sequences
     not see in other reads.


    """
    fout.write(report)

proc summarize(filterLog: string, fn: string) =
    var readsToFilterSum = initTable[string, int]()
    var fstream = newFileStream(fn, fmRead)
    fstream.unpack(readsToFilterSum)
    summarize(filterLog, readsToFilterSum)

proc gapInCoverage(ovls: seq[overlap], minCov: int, minIdt: float): bool =
    ##caculates the coverage in a linear pass. If the start or end < minCov there
##is a gap. The first and last position are skipped
    type PosAndTag = tuple
        pos: int
        tag: bool #starts are true ends are false

    let slop = 50
    var positions = newSeq[PosAndTag]()
    var ep = ovls[0].l1

    for i in ovls:
        if i.idt < minIdt:
            continue
        positions.add( (i.start1 + slop, true))
        positions.add( (i.end1 - slop, false))

    positions.sort()

    var posCoverage = initTable[int, int]()

    var cov = 0
    for i in positions:
        if i.tag == true:
            inc(cov)
        else:
            dec(cov)
        discard posCoverage.hasKeyOrPut(i.pos, 0)
        posCoverage[i.pos] = cov

    #need a little buffer so counts can climb, 500bp from start and end
    for k, v in posCoverage:
        if v < minCov and ((ep - k) > 500) and (k > 500):
            return true
    return false


proc stage1Filter*(overlaps: seq[overlap],
maxDiff: int,
maxOvlp: int,
minOvlp: int,
minLen: int,
minDepth: int,
gapFilt: bool,
minIdt: float,
readsToFilter: var Table[string, int]) =

    var fivePrimeCount, threePrimeCount: int = 0
    var ridA = overlaps[0].ridA
    var containedBreads = initHashSet[string]()
    var aReadPass = true

    if gapFilt:
        if gapInCoverage(overlaps, minDepth, minIdt):
            discard readsToFilter.hasKeyOrPut(ridA, 0)
            readsToFilter[ridA] = readsToFilter[ridA] or GREAD

    for i in overlaps:
        if i.idt < minIdt:
            continue
        if(i.l1 < minLen):
            discard readsToFilter.hasKeyOrPut(i.ridA, 0)
            readsToFilter[i.ridA] = readsToFilter[i.ridA] or SREAD
        if(i.l2 < minLen):
            discard readsToFilter.hasKeyOrPut(i.ridB, 0)
            readsToFilter[i.ridB] = readsToFilter[i.ridB] or SREAD
        if (i.l1 < minLen) or (i.l2 < minLen):
            continue
        if i.tag == "contains":
            containedBreads.incl(i.ridB)
        if i.start1 == 0:
            inc(fivePrimeCount)
        if i.end1 == i.l1:
            inc(threePrimeCount)

    if abs(fivePrimeCount - threePrimeCount) > maxDiff:
        discard readsToFilter.hasKeyOrPut(ridA, 0)
        readsToFilter[ridA] = readsToFilter[ridA] or BREAD
        aReadPass = false
    if fivePrimeCount > maxOvlp:
        discard readsToFilter.hasKeyOrPut(ridA, 0)
        readsToFilter[ridA] = readsToFilter[ridA] or HREAD
        aReadPass = false
    if threePrimeCount > maxOvlp:
        discard readsToFilter.hasKeyOrPut(ridA, 0)
        readsToFilter[ridA] = readsToFilter[ridA] or HREAD
        aReadPass = false
    if fivePrimeCount < minOvlp:
        discard readsToFilter.hasKeyOrPut(ridA, 0)
        readsToFilter[ridA] = readsToFilter[ridA] or LREAD
        aReadPass = false
    if threePrimeCount < minOvlp:
        discard readsToFilter.hasKeyOrPut(ridA, 0)
        readsToFilter[ridA] = readsToFilter[ridA] or LREAD
        aReadPass = false

    if aReadPass:
        for i in containedBreads:
            discard readsToFilter.hasKeyOrPut(i, 0)
            readsToFilter[i] = readsToFilter[i] or CREAD

proc comp(x, y: ScoredOverlap): int =
    ## custom sort for bestN
    if x.ovlpLen == y.ovlpLen:
        if x.mRange < y.mRange:
            return -1
        else:
            return 1
        if x.ovlpLen < y.ovlpLen:
            return -1
        return 1

proc stage2Filter(overlaps: seq[overlap],
 minIdt: float,
 bestN: int,
 readsToFilter: var Table[string, int]): seq[string] =

    var left, right: seq[ScoredOverlap]
    result = newSeq[string]()

    if overlaps[0].ridA in readsToFilter:
        return

    for i in overlaps:
        if i.idt < minIdt:
            continue
        if i.ridB in readsToFilter:
            continue
        if i.start1 == 0:
            left.add((i.score, i.l2 - (i.end2 - i.start2), i))
        elif i.end1 == i.l1:
            right.add((i.score, i.l2 - (i.end2 - i.start2), i))

    left.sort(comp)
    right.sort(comp)

    for i in 0 .. (left.len()-1):
        result.add($left[i].o)
        if i >= bestN and left[i].mRange > 1000:
            break

    for i in 0 .. (right.len()-1):
        result.add($right[i].o)
        if i >= bestN and right[i].mRange > 1000:
            break

proc mergeBlacklists*(blistFofn: string,
 readsToFilter: var Table[string, int]) =
    let f = open(blistFofn, fmRead)
    defer: f.close()

    var bFile: string

    while f.readLine(bFile):
        var tmp = initTable[string, int]()
        var fstream = newFileStream(bFile, fmRead)
        echo "[INFO] merging blacklist file: {bFile}".fmt
        fstream.unpack(tmp)
        for k, v in tmp:
            if k in readsToFilter:
                readsToFilter[k] = (readsToFilter[k] or v)
            else:
                readsToFilter[k] = v

proc dumpBlacklist*(readsToFilter: var Table[string, int]) =
    for k, v in readsToFilter:
        echo "{k:09} {v}".fmt

type
    Stage1 = ref object
        stream: FileStream  # not yet used
        maxDiff: int
        maxCov: int
        minCov: int
        minLen: int
        minIdt: float
        gapFilt: bool
        minDepth: int
        blacklist: string
    Stage2 = ref object
        stream: FileStream  # not yet used
        minIdt: float
        bestN: int
        blacklistIn: string
        filteredOutput: string
proc doStage1(args: Stage1) =
    var readsToFilter1 = initTable[string, int]()
    for i in abOvl():
        stage1Filter(i, args.maxDiff, args.maxCov, args.minCov, args.minLen, args.minDepth, args.gapFilt,
                args.minIdt, readsToFilter1)
    var fstream = newFileStream(args.blacklist, fmWrite)
    fstream.pack(readsToFilter1)
proc doStage2(args: Stage2) =
    var readsToFilter2 = initTable[string, int]()
    var fstream = newFileStream(args.blacklistIn, fmRead)
    fstream.unpack(readsToFilter2)

    var output = open(args.filteredOutput, fmWrite)
    defer: output.close()

    for i in abOvl():
        let lines = stage2Filter(i, args.minIdt, args.bestN, readsToFilter2)
        for l in lines:
            output.writeLine(l)
proc runStage1*(
 maxDiff: int = 100,
 maxCov: int = 200,
 minCov: int = 2,
 minLen: int = 6000,
 minIdt: float = 95.0,
 gapFilt: bool = false,
 minDepth: int = 2,
 blacklist: string
    ) =
    var args = Stage1(
        maxDiff:maxDiff,
        maxCov:maxCov,
        minCov:minCov,
        minLen:minLen,
        minIdt:minIdt,
        gapFilt:gapFilt,
        minDepth:minDepth,
        blacklist:blacklist)
    doStage1(args)

proc runStage2*(
 minIdt: float = 90.0,
 bestN: int = 10,
 blacklistIn: string,
 filteredOutput: string
    ) =
    var args = Stage2(
        minIdt:minIdt,
        bestN:bestN,
        filteredOutput:filteredOutput,
        blacklistIn:blacklistIn)
    doStage2(args)

proc runMergeBlacklists*(blistFofn: string, outFn: string) =
    var readsToFilter = initTable[string, int]()

    mergeBlacklists(blistFofn, readsToFilter)

    var fstream = newFileStream(outFn, fmWrite)
    fstream.pack(readsToFilter)

proc runDumpBlacklist*(blacklist: string) =
    var readsToFilter = initTable[string, int]()
    var fstream = newFileStream(blacklist, fmRead)
    fstream.unpack(readsToFilter)

    dumpBlacklist(readsToFilter)

proc falconRunner*(db: string,
 lasJson: string,
 idtStage1: float = 90.0,
 idtStage2: float = 90.0,
 minLen: int = 6000,
 minCov: int = 2,
 maxCov: int = 200,
 maxDiff: int = 100,
 bestN: int = 10,
 minDepth: int = 2,
 gapFilt: bool = false,
 nProc: int = 24,
 filterLog: string,
 outputFn: string) =
    ##Runs the multi-stage m4 overlap filter for falcon. In stage one, reads that
##trigger a filter are marked including containment, gaps in coverage along the
##A-read, and repeat reads. In stage two the filters are applied and the N-best
##overlaps are kept for the 5prime and 3prime of each read.

    #LasJson is generated by the falcon pipeline
    let js = parseFile(lasJson)
    var stage1 = newSeq[string]()
    var stage2 = newSeq[string]()
    var ovls = newSeq[string]()

    var gapFiltOpts = ""
    if gapFilt:
        gapFiltOpts = "--minDepth {minDepth} --gapFilt".fmt

    var counter = 0
    var msgpckFofnS1 = open("blacklist_msgpck.stage1.fofn", fmWrite)

    for las in js:
        inc(counter)
        stage1.add(


                                     "time LA4Falcon -mo {db} {las.getStr()} | falconc m4filt-stage1 --minIdt {idtStage1} {gapFiltOpts} --minLen {minLen} --minCov {minCov} --maxCov {maxCov} --maxDiff {maxDiff} --blacklist {counter}.stage1.tmp.msgpck".fmt)
        stage2.add(


                                     "time LA4Falcon -mo {db} {las.getStr()} | falconc m4filt-stage2 --minIdt {idtStage2} --bestN {bestN}  --blacklistIn merged_blacklist.stage1.msgpck --filteredOutput {counter}.tmp.ovl".fmt)
        ovls.add("{counter}.tmp.ovl".fmt)
        msgpckFofnS1.writeLine("{counter}.stage1.tmp.msgpck".fmt)
    close(msgpckFofnS1)

    let stage1rv = execProcesses(stage1, options = {poEchoCmd}, n = nProc)
    if stage1rv != 0:
        quit "[FATAL] stage one had a failure!"

    let merge1 = execCmd(



                                 "falconc m4filt-merge-blacklist -b blacklist_msgpck.stage1.fofn -o merged_blacklist.stage1.msgpck")
    if merge1 != 0:
        quit "[FATAL] merge one failure!"

    let stage2rv = execProcesses(stage2, options = {poEchoCmd}, n = nProc)
    if stage2rv != 0:
        quit "[FATAL] stage two had a failure!"

    var outFile = open(outputFn, fmWrite)
    defer: outFile.close()

    for i in ovls:
        let ov = open(i, fmRead)
        defer: ov.close()
        for line in ov.lines:
            outFile.writeLine(line)
    outFile.writeLine("---")

    summarize(filterLog, "merged_blacklist.stage1.msgpck")

proc ipaRunner*(ovlsFofn: string,
 idtStage1: float = 90.0,
 idtStage2: float = 90.0,
 minLen: int = 6000,
 minCov: int = 2,
 maxCov: int = 200,
 maxDiff: int = 100,
 bestN: int = 10,
 minDepth: int = 2,
 gapFilt: bool = false,
 nProc: int = 24,
 filterLog: string,
 outputFn: string) =
    ##Runs the multi-stage m4 overlap filter for IPA. In stage one, reads that
    ##trigger a filter are marked including containment, gaps in coverage along
    ##the A-read, and repeat reads. In stage two the filters are applied and the
    ##N-best overlaps are kept for the 5prime and 3prime of each read.

    #fofn of m4 files generated in the IPA WDL
    let fh = open(ovlsFofn, fmRead)
    var m4s = newSeq[string]()
    for l in fh.lines:
        if isEmptyFile(l):
            continue
        m4s.add(l)

    var stage1 = newSeq[string]()
    var stage2 = newSeq[string]()
    var ovls = newSeq[string]()

    var gapFiltOpts = ""
    if gapFilt:
        gapFiltOpts = "--minDepth {minDepth} --gapFilt".fmt
    var counter = 0
    var msgpckFofnS1 = open("blacklist_msgpck.stage1.fofn", fmWrite)
    for ms in m4s:
        inc(counter)
        stage1.add("cat {ms} | falconc m4filt-stage1 --minIdt {idtStage1} {gapFiltOpts} --minLen {minLen} --minCov {minCov} --maxCov {maxCov} --maxDiff {maxDiff} --blacklist {counter}.stage1.tmp.msgpck".fmt)
        stage2.add("cat {ms} | falconc m4filt-stage2 --minIdt {idtStage2} --bestN {bestN}  --blacklistIn merged_blacklist.stage1.msgpck --filteredOutput {counter}.tmp.ovl".fmt)
        ovls.add("{counter}.tmp.ovl".fmt)
        msgpckFofnS1.writeLine("{counter}.stage1.tmp.msgpck".fmt)
    close(msgpckFofnS1)

    let stage1rv = execProcesses(stage1, options = {poEchoCmd}, n = nProc)
    if stage1rv != 0:
        quit "[FATAL] stage one had a failure!"

    let merge1 = execCmd("falconc m4filt-merge-blacklist -b blacklist_msgpck.stage1.fofn -o merged_blacklist.stage1.msgpck")
    if merge1 != 0:
        quit "[FATAL] merge one failure!"

    let stage2rv = execProcesses(stage2, options = {poEchoCmd}, n = nProc)
    if stage2rv != 0:
        quit "[FATAL] stage two had a failure!"

    var outFile = open(outputFn, fmWrite)
    defer: outFile.close()
    for i in ovls:
        let ov = open(i, fmRead)
        defer: ov.close()
        for line in ov.lines:
            outFile.writeLine(line)
    outFile.writeLine("---")

    summarize(filterLog, "merged_blacklist.stage1.msgpck")


proc main() =
    echo "main"

when isMainModule:
    main()
