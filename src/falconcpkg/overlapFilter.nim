import threadpool_simple
import sets
import tables
from strutils import split, parseInt, parseFloat, parseBool
#from sequtils import keepIf
from strformat import fmt
from ./util import isEmptyFile, log
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

proc parseOvl*(s: string): overlap =
    var ovl: overlap
    let ld = s.split(" ")
    if 13 != ld.len():
        let msg = "Error parsing ovl (split={ld.len()}): '{s}'".fmt
        log(msg)
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



iterator abOvl*(sin: Stream): seq[overlap] =
    var
        readA: string
        ov: overlap
        #The buffer len of 500 should be adjusted if more fields are added.
        buff = newString(500)
        ovls = newSeq[overlap]()
    discard readLine(sin, buff)
    try:
        ov = parseOvl(buff)
        readA = ov.ridA
        ovls.add(ov)
        while(readLine(sin, buff)):
            ov = parseOvl(buff)
            if ov.ridA != readA:
                yield ovls
            if ov.ridA != readA:
                readA = ov.ridA
                ovls = @[]
            ovls.add(ov)
        yield ovls
    except AssertionError:
        log("buff='{buff}'".fmt)
        raise



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
    defer: fout.close()
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
    defer: fstream.close()
    fstream.unpack(readsToFilterSum)
    summarize(filterLog, readsToFilterSum)

proc gapInCoverage*(ovls: seq[overlap], minDepth: int, minIdt: float): bool =
    ##Calculates the coverage in a linear pass. If the start or end < minDepth there
    ##is a gap. The first and last position are skipped
    type PosAndTag = tuple
        pos: int #read a position
        tag: bool #starts are true ends are false
        cli: bool #if start == 0 or end == len; is it a clipped read
        loc: bool #is this a local align
    type Info = tuple
        cov: int #coverage
        cco: int #clean cov
        cli: int #clip


    var ep = ovls[0].l1
    var an = ovls[0].ridA

    var positions = newSeq[PosAndTag]()
    # load start/end into a tuple for linear depth caculation
    for i in ovls:
        if i.idt < 95.0:
            continue
        var clipped = (i.start2 != 0) and (i.end2 != i.l2)
        var loc = (i.tag == "u")
        positions.add( (i.start1, true, clipped, loc))
        positions.add( (i.end1, false,  clipped, loc))

    positions.sort()

    var runningClip, runningCov, runningClean = 0
    var posInfo = initOrderedTable[int, Info]()
    # turn running start/end into a depth at each start/end
    for i in 0 .. (positions.len() - 1):
        if positions[i].tag:
            inc(runningCov)
            if not positions[i].loc:
                inc(runningClean)
            if positions[i].cli:
                inc(runningClip)
        else:
            dec(runningCov)
            dec(runningClip)
            if not positions[i].loc:
                dec(runningClean)
        posInfo[positions[i].pos] = (runningCov, runningClean, max(0,
                runningClip))


    var hasCovDip = false
    var hasZeroCovInDip = false
    var clipHigh = false
    var lastOkCov = 0
    var lastOkPos = 0
    var lastZero = -1
    var lastPos   = 0
    var lastCleanCov = 0
    var count = 0
    var averageCoverage: float  = 0

    var clipCans = newSeq[int]()
    var depthCans = newSeq[int]()
    for i, j in posInfo:
        #calculating depth sum over the lengths
        if count != (posInfo.len - 1) and count != 0:
            averageCoverage += float((i - lastPos) * lastCleanCov)
        lastPos = i
        lastCleanCov = j.cco
        # clipping coverage is almost as high/higher than clean coverage
        if (j.cco - j.cli) < 5 and (j.cli > 2):
            clipHigh = true
            clipCans.add(i)
        if j.cco == 0:
            lastZero = i
        if j.cco >= minDepth:
            if (count - lastOkCov) > 1 and (lastOkCov != 0):
                hasCovDip = true
                if lastOkPos < lastZero and i > lastZero:
                 hasZeroCovInDip = true
                depthCans.add(i)
                depthCans.add(lastOkPos)
            lastOkCov = count
            lastOkPos = i
        inc(count)
    # Please leave these print statements, they are critical for troubleshooting
    #    echo "{an} {i} {j.cov} allCov".fmt
    #    echo "{an} {i} {j.cco} cleanCov".fmt
    #    echo "{an} {i} {j.cli} clipCov".fmt
    averageCoverage = averageCoverage/float(ep)
    #echo "clipHigh:{clipHigh} hasCovDip:{hasCovDip} hasZeroDepthInDip:{hasZeroCovInDip} avgCov:{averageCoverage}".fmt
    if averageCoverage < 5:
        return false
    # Requested by Ivan
    if hasCovDip and hasZeroCovInDip:
        return true
    if clipHigh and hasCovDip:
        for dp in depthCans:
            for cp in clipCans:
                #stderr.writeLine("POSS: {dp} {cp}".fmt)
                if abs(dp - cp) < 500:
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
            #stderr.writeLine("YES {ridA}".fmt)

    for i in overlaps:
        if i.tag == "u":
            continue
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
        if i.tag == "u":
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
        defer: fstream.close()
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
        icmd: string
        sin: Stream
        maxDiff: int
        maxCov: int
        minCov: int
        minLen: int
        minIdt: float
        gapFilt: bool
        minDepth: int
        blacklist: string
    Stage2 = ref object
        icmd: string
        sin: Stream
        minIdt: float
        bestN: int
        blacklistIn: string
        filteredOutput: string
proc doStage1(args: Stage1) =
    var readsToFilter1 = initTable[string, int]()
    for i in abOvl(args.sin):
        stage1Filter(i, args.maxDiff, args.maxCov, args.minCov, args.minLen,
                args.minDepth, args.gapFilt,
                args.minIdt, readsToFilter1)
    var fstream = newFileStream(args.blacklist, fmWrite)
    defer: fstream.close()
    fstream.pack(readsToFilter1)
proc doStage2(args: Stage2) =
    var readsToFilter2 = initTable[string, int]()
    var fstream = newFileStream(args.blacklistIn, fmRead)
    defer: fstream.close()
    fstream.unpack(readsToFilter2)

    var output = open(args.filteredOutput, fmWrite)
    defer: output.close()

    for i in abOvl(args.sin):
        let lines = stage2Filter(i, args.minIdt, args.bestN, readsToFilter2)
        for l in lines:
            output.writeLine(l)
proc startStage1(args: Stage1) =
    log("startStage1: ", args.icmd)
    var p = startProcess(command = args.icmd, options = {poEvalCommand,
            })
            #poStdErrToStdOut})
    var argsx = args
    argsx.sin = p.outputStream
    doStage1(argsx)
    p.close()
proc startStage2(args: Stage2) =
    log("startStage2: ", args.icmd)
    var p = startProcess(command = args.icmd, options = {poEvalCommand,
            })
            #poStdErrToStdOut})
    var argsx = args
    argsx.sin = p.outputStream
    doStage2(argsx)
    p.close()
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
        maxDiff: maxDiff,
        maxCov: maxCov,
        minCov: minCov,
        minLen: minLen,
        minIdt: minIdt,
        gapFilt: gapFilt,
        minDepth: minDepth,
        blacklist: blacklist)
    args.sin = newFileStream(stdin)
    doStage1(args)

proc runStage2*(
 minIdt: float = 90.0,
 bestN: int = 10,
 blacklistIn: string,
 filteredOutput: string
    ) =
    var args = Stage2(
        minIdt: minIdt,
        bestN: bestN,
        filteredOutput: filteredOutput,
        blacklistIn: blacklistIn)
    args.sin = newFileStream(stdin)
    doStage2(args)

proc runMergeBlacklists*(blistFofn: string, outFn: string) =
    var readsToFilter = initTable[string, int]()

    mergeBlacklists(blistFofn, readsToFilter)

    var fstream = newFileStream(outFn, fmWrite)
    defer: fstream.close()
    fstream.pack(readsToFilter)

proc runDumpBlacklist*(blacklist: string) =
    var readsToFilter = initTable[string, int]()
    var fstream = newFileStream(blacklist, fmRead)
    defer: fstream.close()
    fstream.unpack(readsToFilter)

    dumpBlacklist(readsToFilter)

type
    M4filtOptions = object
        idtStage1: float
        idtStage2: float
        minLen: int
        minCov: int
        maxCov: int
        maxDiff: int
        bestN: int
        minDepth: int
        gapFilt: bool
        keepIntermediates: bool
        filterLogFn: string
        outFn: string

proc m4filt(icmds: seq[string],
    opts: M4filtOptions,
    threadpool: threadpool_simple.ThreadPool) =

    var stage1 = newSeq[Stage1]()
    var stage2 = newSeq[Stage2]()
    var ovls = newSeq[string]()

    let minDepthGapFilt =
        if opts.gapFilt:
            opts.minDepth
        else:
            0

    var counter = 0
    let blacklist_msgpack_fn = "blacklist_fofn.stage1.msgpck"
    var msgpackFofnS1 = open(blacklist_msgpack_fn, fmWrite)

    let merged_blacklist_msgpack_fn = "merged_blacklist.stage1.msgpck"

    var intermediateFns: seq[string]
    intermediateFns.add(blacklist_msgpack_fn)
    intermediateFns.add(merged_blacklist_msgpack_fn)

    for icmd in icmds:
        inc(counter)

        let blacklist_fn = "{counter}.stage1.tmp.msgpck".fmt
        intermediateFns.add(blacklist_fn)
        let args1 = Stage1(
            icmd: icmd,
            maxDiff: opts.maxDiff,
            maxCov: opts.maxCov,
            minCov: opts.minCov,
            minLen: opts.minLen,
            minIdt: opts.idtStage1,
            gapFilt: opts.gapFilt,
            minDepth: minDepthGapFilt,
            blacklist: blacklist_fn)
        stage1.add(args1)

        # We will generate merged_blacklist_msgpack btw stage 1 and 2.
        let ovlFn = "{counter}.tmp.ovl".fmt
        let args2 = Stage2(
            icmd: icmd,
            minIdt: opts.idtStage2,
            bestN: opts.bestN,
            filteredOutput: ovlFn,
            blacklistIn: merged_blacklist_msgpack_fn)
        stage2.add(args2)

        ovls.add(ovlFn)
        msgpackFofnS1.writeLine(blacklist_fn)
    msgpackFofnS1.close()

    for x in stage1:
        threadpool.spawn startStage1(x)
        #startStage1(x)
    threadpool.sync()

    runMergeBlacklists(blacklist_msgpack_fn, merged_blacklist_msgpack_fn)
    summarize(opts.filterLogFn, merged_blacklist_msgpack_fn)

    for x in stage2:
        threadpool.spawn startStage2(x)
        #startStage2(x)
    threadpool.sync()

    var outFile = open(opts.outFn, fmWrite)

    for i in ovls:
        let ov = open(i, fmRead)
        defer: ov.close()
        for line in ov.lines:
            outFile.writeLine(line)
    outFile.writeLine("---")
    outFile.close()

    if not opts.keepIntermediates:
        util.removeFiles(ovls)
        util.removeFiles(intermediateFns)

proc m4filtRunner*(
 idtStage1: float = 90.0,
 idtStage2: float = 90.0,
 minLen: int = 6000,
 minCov: int = 2,
 maxCov: int = 200,
 maxDiff: int = 100,
 bestN: int = 10,
 minDepth: int = 2,
 gapFilt: bool = false,
 keepIntermediates: bool = false,
 nProc: int = 24,
 filterLogFn: string,
 outFn: string) =
    ## Run the multi-stage m4 overlap filter.
    ## Read the m4 filenames from stdin.
    ## In stage one, reads that
    ## trigger a filter are marked including containment, gaps in coverage along the
    ## A-read, and repeat reads.
    ## In stage two the filters are applied and the N-best
    ## overlaps are kept for the 5prime and 3prime of each read.

    var icmds: seq[string]
    for line in lines(stdin):
        let icmd = "cat {line}".fmt
        icmds.add(icmd)

    let opts = M4filtOptions(
        idtStage1: idtStage1,
        idtStage2: idtStage2,
        minLen: minLen,
        minCov: minCov,
        maxCov: maxCov,
        maxDiff: maxDiff,
        bestN: bestN,
        minDepth: minDepth,
        gapFilt: gapFilt,
        keepIntermediates: keepIntermediates,
        filterLogFn: filterLogFn,
        outFn: outFn)

    let threadpool = threadpool_simple.newThreadPool(nProc)

    m4filt(icmds, opts, threadpool)

proc falconRunner*(db: string,
 lasJsonFn: string,
 idtStage1: float = 90.0,
 idtStage2: float = 90.0,
 minLen: int = 6000,
 minCov: int = 2,
 maxCov: int = 200,
 maxDiff: int = 100,
 bestN: int = 10,
 minDepth: int = 2,
 gapFilt: bool = false,
 keepIntermediates: bool = false,
 nProc: int = 24,
 filterLogFn: string,
 outFn: string) =
    ##Runs the multi-stage m4 overlap filter for falcon. In stage one, reads that
    ##trigger a filter are marked including containment, gaps in coverage along the
    ##A-read, and repeat reads. In stage two the filters are applied and the N-best
    ##overlaps are kept for the 5prime and 3prime of each read.

    #LasJson is generated by the falcon pipeline
    let js = parseFile(lasJsonFn)
    var icmds: seq[string]
    for las in js:
        # Remember to redirect stderr, since we will not read it.
        # Otherwise, beyond several thousand lines, the stderr buffer
        # will fill up and block writing to stdout.
        let icmd = "LA4Falcon -mo {db} {las.getStr()} 2> /dev/null".fmt
        icmds.add(icmd)
    let opts = M4filtOptions(
        idtStage1: idtStage1,
        idtStage2: idtStage2,
        minLen: minLen,
        minCov: minCov,
        maxCov: maxCov,
        maxDiff: maxDiff,
        bestN: bestN,
        minDepth: minDepth,
        gapFilt: gapFilt,
        keepIntermediates: keepIntermediates,
        filterLogFn: filterLogFn,
        outFn: outFn)

    let threadpool = threadpool_simple.newThreadPool(nProc)

    m4filt(icmds, opts, threadpool)

proc ipaRunner*(ovlsFofnFn: string,
 idtStage1: float = 90.0,
 idtStage2: float = 90.0,
 minLen: int = 6000,
 minCov: int = 2,
 maxCov: int = 200,
 maxDiff: int = 100,
 bestN: int = 10,
 minDepth: int = 2,
 gapFilt: bool = false,
 keepIntermediates: bool = false,
 nProc: int = 24,
 filterLogFn: string,
 outputFn: string) =
    ##Runs the multi-stage m4 overlap filter for IPA. In stage one, reads that
    ##trigger a filter are marked including containment, gaps in coverage along
    ##the A-read, and repeat reads. In stage two the filters are applied and the
    ##N-best overlaps are kept for the 5prime and 3prime of each read.

    #fofn of m4 files generated in the IPA WDL
    let fh = open(ovlsFofnFn, fmRead)
    defer: fh.close()
    var m4s = newSeq[string]()
    for l in fh.lines:
        if isEmptyFile(l):
            continue
        m4s.add(l)

    var icmds: seq[string]
    for ms in m4s:
        let icmd = "/bin/cat {ms}".fmt
        icmds.add(icmd)

    let opts = M4filtOptions(
        idtStage1: idtStage1,
        idtStage2: idtStage2,
        minLen: minLen,
        minCov: minCov,
        maxCov: maxCov,
        maxDiff: maxDiff,
        bestN: bestN,
        minDepth: minDepth,
        gapFilt: gapFilt,
        keepIntermediates: keepIntermediates,
        filterLogFn: filterLogFn,
        outFn: outputFn)

    let threadpool = threadpool_simple.newThreadPool(nProc)

    m4filt(icmds, opts, threadpool)


proc main() =
    echo "main"

when isMainModule:
    main()
