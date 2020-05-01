# vim: sw=4 ts=4 sts=4 tw=0 et:
import times
import threadpool_simple
import sets
import tables
#from sequtils import keepIf
from strformat import fmt
#from os import getFileInfo
from system/io import getFileSize
from ./util import isEmptyFile, log
from ./overlapParser import Overlap, parseOverlap
import ./overlapParser as op
#, getNextPile, parseM4, toString, index, M4Index
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
 13. tag (optional)
 14. extra column(s)? (optional)
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
    coverageInfo = tuple[gap: bool, lowCov: bool, highCov: bool,
            balance: bool]

    ScoredOverlap = tuple
        ovlpLen: int
        mRange: int
        o: Overlap

proc lowIdt*(o: Overlap, idt: float): bool =
    ##Check if percent identity is below threshold
    #This is a overlap filter
    if o.idt < idt:
        return true
    return false

proc checkFractionOverlap*(o: Overlap, lowOverlap: float,
        highOverlap: float): bool =
    ##Check for percentage of two reads that overlap
    #This is an overlap filter - new as of June 2019
    let perA = (float(o.Aend - o.Astart) / (float(o.Alen))) * 100
    let perB = (float(o.Bend - o.Bstart) / (float(o.Blen))) * 100
    if (perA < lowOverlap) or (perA > highOverlap) or (perB < lowOverlap) or (
            perB > highOverlap):
        return true
    return false

proc small*(o: Overlap, l: int): bool =
    ##Tests if a read is too small
    #This is a read filter
    if o.Alen < l:
        return true
    else:
        return false

proc contained*(o: Overlap): bool =
    ##Tests if A-read is contained within B
    #This is a read filter
    if o.Astart > 0:
        return false
    if o.Aend != o.Alen:
        return false
    else:
        return true

proc missingTerminus*(o: Overlap): bool =
    ##Test that overlap reaches the end of both reads
    #This is an overlap filter
    if o.Astart != 0 and o.Aend != o.Alen:
        return true
    if o.Bstart != 0 and o.Bend != o.Blen:
        return true
    else:
        return false



proc summarize(filterLog: string, readsToFilter: var tables.Table[string, int]) =
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
    var readsToFilterSum = tables.initTable[string, int]()
    var fstream = newFileStream(fn, fmRead)
    defer: fstream.close()
    fstream.unpack(readsToFilterSum)
    summarize(filterLog, readsToFilterSum)

proc gapInCoverage*(ovls: seq[Overlap], minDepth: int, minIdt: float): bool =
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


    var ep = ovls[0].Alen
    var an = ovls[0].Aname

    var positions = newSeq[PosAndTag]()
    # load start/end into a tuple for linear depth caculation
    for i in ovls:
        if i.idt < 95.0:
            continue
        var clipped = (i.Bstart != 0) and (i.Bend != i.Blen)
        var loc = (i.tag == "u")
        positions.add( (i.Astart, true, clipped, loc))
        positions.add( (i.Aend, false, clipped, loc))

    positions.sort()

    var runningClip, runningCov, runningClean = 0
    var posInfo = tables.initOrderedTable[int, Info]()
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
    var lastPos = 0
    var lastCleanCov = 0
    var count = 0
    var averageCoverage: float = 0

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

proc stage1Filter*(overlaps: seq[Overlap],
 maxDiff: int,
 maxOvlp: int,
 minOvlp: int,
 minLen: int,
 minDepth: int,
 gapFilt: bool,
 minIdt: float,
 readsToFilter: var tables.Table[string, int]) =
    if 0 == len(overlaps):
        return
    var fivePrimeCount, threePrimeCount: int = 0
    var ridA = overlaps[0].Aname
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
        if(i.Alen < minLen):
            discard readsToFilter.hasKeyOrPut(i.Aname, 0)
            readsToFilter[i.Aname] = readsToFilter[i.Aname] or SREAD
        if(i.Blen < minLen):
            discard readsToFilter.hasKeyOrPut(i.Bname, 0)
            readsToFilter[i.Bname] = readsToFilter[i.Bname] or SREAD
        if (i.Alen < minLen) or (i.Blen < minLen):
            continue
        if i.tag == "contains":
            containedBreads.incl(i.Bname)
        if i.Astart == 0:
            inc(fivePrimeCount)
        if i.Aend == i.Alen:
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

proc CheckOverhangLen*(ovl: Overlap, minOverhang: int): bool =
    ## This function should filter only dovetail overlaps which
    ## do not extend beyond minOverhang bases over left or
    ## right ends.
    ## Any non-dovetail overlap should return true, because this function
    ## doesn't filter those.

    let a_left = ovl.Astart
    let a_right = ovl.Alen - ovl.Aend
    var b_left = ovl.Bstart
    var b_right = ovl.Blen - ovl.Bend
    if ovl.Brev:
        swap(b_left, b_right)
    let a_hang = a_left - b_left
    let b_hang = b_right - a_right

    let allowedDovetailDist: int = 0

    # Unknown.
    var ovlType: string = "u"

    if a_left <= allowedDovetailDist and a_right <= allowedDovetailDist:
        # A is contained in B.
        ovlType = "c"
    elif b_left <= allowedDovetailDist and b_right <= allowedDovetailDist:
        # B is contained in A.
        ovlType = "C"
    elif a_left <= allowedDovetailDist and b_right <= allowedDovetailDist:
        # 5' overlap.
        ovlType = "5"
    elif a_right <= allowedDovetailDist and b_left <= allowedDovetailDist:
        # 3' overlap.
        ovlType = "3"

    if ovlType != "3" and ovlType != "5":
        return true

    # Legend:
    # http://wgs-assembler.sourceforge.net/wiki/index.php/Overlaps
    #   - 3' overlap: a_hang >= 0 && b_hang >= 0
    #   - 5' overlap: a_hang <= 0 && b_hang <= 0
    #   - B contained in A: a_hang >= 0 && b_hang <= 0
    #   - A contained in B: a_hang <= 0 && b_hang >= 0
    if a_hang >= 0 and b_hang >= 0 and a_hang < minOverhang and b_hang < minOverhang:
        return false
    elif a_hang <= 0 and b_hang <= 0 and a_hang > (-minOverhang) and b_hang > (-minOverhang):
        return false

    return true

proc stage2Filter(overlaps: seq[Overlap],
 minIdt: float,
 bestN: int,
 minOverhang: int,
 readsToFilter: var tables.Table[string, int]): seq[string] =

    var left, right: seq[ScoredOverlap]
    result = newSeq[string]()

    if 0 == len(overlaps):
        return
    assert len(overlaps) > 0
    if overlaps[0].Aname in readsToFilter:
        return

    for i in overlaps:
        if i.idt < minIdt:
            continue
        if i.Bname in readsToFilter:
            continue
        if i.tag == "u":
            continue
        if minOverhang > 0 and CheckOverhangLen(i, minOverhang) == false:
            continue

        if i.Astart == 0:
            left.add((i.score, i.Blen - (i.Bend - i.Bstart), i))
        elif i.Aend == i.Alen:
            right.add((i.score, i.Blen - (i.Bend - i.Bstart), i))

    left.sort(comp)
    right.sort(comp)

    for i in 0 .. (left.len()-1):
        result.add(toString(left[i].o))
        if i >= bestN and left[i].mRange > 1000:
            break

    for i in 0 .. (right.len()-1):
        result.add(toString(right[i].o))
        if i >= bestN and right[i].mRange > 1000:
            break

proc mergeBlacklists*(blistFofn: string,
 readsToFilter: var tables.Table[string, int]) =
    let f = open(blistFofn, fmRead)
    defer: f.close()

    var bFile: string

    while f.readLine(bFile):
        var tmp = tables.initTable[string, int]()
        var fstream = newFileStream(bFile, fmRead)
        defer: fstream.close()
        log("merging blacklist file: {bFile}".fmt)
        fstream.unpack(tmp)
        for k, v in tmp:
            if k in readsToFilter:
                readsToFilter[k] = (readsToFilter[k] or v)
            else:
                readsToFilter[k] = v

proc dumpBlacklist*(readsToFilter: var tables.Table[string, int]) =
    for k, v in readsToFilter:
        echo "{k:09} {v}".fmt

type
    Stage1 = ref object
        maxDiff: int
        maxCov: int
        minCov: int
        minLen: int
        minIdt: float
        gapFilt: bool
        minDepth: int
        blacklist: string
        # Used only san M4Index:
        icmd: string
        sin: Stream
        # Used only with the M4Index:
        index: M4Index
        m4Fn: string
    Stage2 = ref object
        minIdt: float
        bestN: int
        minOverhang: int
        blacklistIn: string
        filteredOutput: string
        # Used only san M4Index:
        icmd: string
        sin: Stream
        # Used only with the M4Index:
        index: M4Index
        m4Fn: string
proc doStage1(args: Stage1) =
    var readsToFilter1 = tables.initTable[string, int]()
    for i in op.getNextPile(args.sin):
        stage1Filter(i, args.maxDiff, args.maxCov, args.minCov, args.minLen,
                args.minDepth, args.gapFilt,
                args.minIdt, readsToFilter1)
    var fstream = newFileStream(args.blacklist, fmWrite)
    defer: fstream.close()
    fstream.pack(readsToFilter1)
proc doStage1Indexed(args: Stage1) =
    # Here, we actually use the index in the Stage1 obj.
    var readsToFilter1 = tables.initTable[string, int]()
    var sin = streams.newFileStream(args.m4Fn, fmRead)
    defer: sin.close()
    for i in getNextPile(sin, args.index):
        stage1Filter(i, args.maxDiff, args.maxCov, args.minCov, args.minLen,
                args.minDepth, args.gapFilt,
                args.minIdt, readsToFilter1)
    var fstream = newFileStream(args.blacklist, fmWrite)
    defer: fstream.close()
    fstream.pack(readsToFilter1)
proc doStage2(args: Stage2) =
    var readsToFilter2 = tables.initTable[string, int]()
    var fstream = newFileStream(args.blacklistIn, fmRead)
    defer: fstream.close()
    fstream.unpack(readsToFilter2)

    var output = open(args.filteredOutput, fmWrite)
    defer: output.close()

    for i in op.getNextPile(args.sin):
        let lines = stage2Filter(i, args.minIdt, args.bestN, args.minOverhang, readsToFilter2)
        for l in lines:
            output.writeLine(l)
proc doStage2Indexed(args: Stage2) =
    # Here, we actually use the index in the Stage2 obj.
    var readsToFilter2 = tables.initTable[string, int]()
    var fstream = newFileStream(args.blacklistIn, fmRead)
    defer: fstream.close()
    fstream.unpack(readsToFilter2)

    var output = open(args.filteredOutput, fmWrite)
    defer: output.close()

    var sin = streams.newFileStream(args.m4Fn, fmRead)
    defer: sin.close()
    for i in getNextPile(sin, args.index):
        let lines = stage2Filter(i, args.minIdt, args.bestN, args.minOverhang, readsToFilter2)
        for l in lines:
            output.writeLine(l)
proc startStage1(args: Stage1) =
    log("startStage1: ", args.icmd)
    var p = osproc.startProcess(command = args.icmd, options = {poEvalCommand,
        })
        #poStdErrToStdOut})
    if osproc.peekExitCode(p) > 0:
        let msg = "Immediate failure in stage1 startProcess('" & args.icmd & "')"
        util.raiseEx(msg)
    var argsx = args
    argsx.sin = p.outputStream
    try:
        doStage1(argsx)
    finally:
        osproc.close(p)
        if osproc.peekExitCode(p) > 0:
            let msg = "Failure in stage2 startProcess('" & args.icmd & "')"
            util.raiseEx(msg)
proc startStage2(args: Stage2) =
    log("startStage2: ", args.icmd)
    var p = osproc.startProcess(command = args.icmd, options = {poEvalCommand,
        })
        #poStdErrToStdOut})
    if osproc.peekExitCode(p) > 0:
        let msg = "Immediate failure in stage2 startProcess('" & args.icmd & "')"
        util.raiseEx(msg)
    var argsx = args
    argsx.sin = p.outputStream
    try:
        doStage2(argsx)
    finally:
        osproc.close(p)
    if osproc.peekExitCode(p) > 0:
        let msg = "Failure in stage2 startProcess('" & args.icmd & "')"
        util.raiseEx(msg)
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
 minOverhang: int = 0,
 blacklistIn: string,
 filteredOutput: string
    ) =
    var args = Stage2(
        minIdt: minIdt,
        bestN: bestN,
        minOverhang: minOverhang,
        filteredOutput: filteredOutput,
        blacklistIn: blacklistIn)
    args.sin = newFileStream(stdin)
    doStage2(args)

proc runMergeBlacklists*(blistFofn: string, outFn: string) =
    var readsToFilter = tables.initTable[string, int]()

    mergeBlacklists(blistFofn, readsToFilter)

    var fstream = newFileStream(outFn, fmWrite)
    defer: fstream.close()
    fstream.pack(readsToFilter)

proc runDumpBlacklist*(blacklist: string) =
    var readsToFilter = tables.initTable[string, int]()
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
        minOverhang: int
        minDepth: int
        gapFilt: bool
        keepIntermediates: bool
        filterLogFn: string
        outFn: string

proc m4filtSingleton(m4Fn: string,
  splits: seq[int],
  index: M4Index,
  opts: M4filtOptions,
  threadpool: threadpool_simple.ThreadPool) =
    # This m4filt operates on a single M4 file.
    # The parallelism is on split portions of it.

    log("Pile-weighted splits:", splits)
    #util.raiseEx("HELLO")
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

    var
        pileBeg = 0
        pileEnd = 0

    for npiles in splits:
        inc(counter)
        pileBeg = pileEnd
        pileEnd = pileBeg + npiles

        let blacklist_fn = "{counter}.stage1.tmp.msgpck".fmt
        intermediateFns.add(blacklist_fn)
        let args1 = Stage1(
            m4Fn: m4Fn,
            index: index[pileBeg ..< pileEnd],
            #icmd: icmd,
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
            m4Fn: m4Fn,
            index: index[pileBeg ..< pileEnd],
            #icmd: icmd,
            minIdt: opts.idtStage2,
            bestN: opts.bestN,
            minOverhang: opts.minOverhang,
            filteredOutput: ovlFn,
            blacklistIn: merged_blacklist_msgpack_fn)
        stage2.add(args2)

        ovls.add(ovlFn)
        msgpackFofnS1.writeLine(blacklist_fn)
    msgpackFofnS1.close()

    let timeStart = times.getTime()
    for x in stage1:
        threadpool.spawn doStage1Indexed(x)
        #doStage1Indexed(x)
    threadpool.sync()

    let timeStage1 = times.getTime()
    log("TIME stage1:", (timeStage1 - timeStart))
    runMergeBlacklists(blacklist_msgpack_fn, merged_blacklist_msgpack_fn)
    summarize(opts.filterLogFn, merged_blacklist_msgpack_fn)

    let timeBlacklist = times.getTime()
    log("TIME blacklist:", (timeBlacklist - timeStage1))
    for x in stage2:
        threadpool.spawn doStage2Indexed(x)
        #doStage2Indexed(x)
    threadpool.sync()
    let timeEnd = times.getTime()

    log("TIME stage2:", (timeEnd - timeBlacklist))

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

proc m4filtPiped(icmds: seq[string],
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
            #index: index,
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

    let timeStart = times.getTime()
    for x in stage1:
        threadpool.spawn startStage1(x)
        #startStage1(x)
    threadpool.sync()

    let timeStage1 = times.getTime()
    log("TIME stage1:", (timeStage1 - timeStart))
    runMergeBlacklists(blacklist_msgpack_fn, merged_blacklist_msgpack_fn)
    summarize(opts.filterLogFn, merged_blacklist_msgpack_fn)

    let timeBlacklist = times.getTime()
    log("TIME blacklist:", (timeBlacklist - timeStage1))
    for x in stage2:
        threadpool.spawn startStage2(x)
        #startStage2(x)
    threadpool.sync()
    let timeEnd = times.getTime()

    log("TIME stage2:", (timeEnd - timeBlacklist))

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

proc m4filtContainedStreams*(
 sin: streams.Stream,
 sout: streams.Stream,
 min_len: int,
 min_idt_pct: float): int =
    # Return the number of lines written.
    let overlaps = op.parseM4(sin)

    # Filter for length and identity
    var good_enough_overlaps: seq[Overlap]

    for ovl in overlaps:
        if ovl.Aname == ovl.Bname: # don't need self-self overlapping
            continue
        if ovl.tag == "none" or ovl.tag == "?":
            continue
        if ovl.idt < min_idt_pct: # only take record with >96% identity as overlapped reads
            continue
        # Only use reads longer than min_len for assembly.
        if ovl.Alen < min_len:
            continue
        if ovl.Blen < min_len:
            continue
        good_enough_overlaps.add(ovl)

    # Find all contained rids.
    var contained_rids = sets.initHashSet[string]()
    for ovl in good_enough_overlaps:
        if ovl.tag == "contained" or ovl.tag == "C":
            contained_rids.incl(ovl.Aname)
        elif ovl.tag == "contains" or ovl.tag == "c":
            contained_rids.incl(ovl.Bname)

    # Drop all overlaps with contained reads.
    var desired_overlaps: seq[Overlap]

    for ovl in good_enough_overlaps:
        if contained_rids.contains(ovl.Aname) or contained_rids.contains(ovl.Bname):
            continue
        desired_overlaps.add(ovl)

    for ovl in desired_overlaps:
        sout.writeLine(toString(ovl))
        result += 1

proc m4filtContained*(
 in_fn: string,
 out_fn: string,
 lfc = false, # IGNORED
 disable_chimer_bridge_removal = false, # IGNORED
 ctg_prefix = "", # IGNORED
 min_len: int = 400,
 min_idt_pct: float = 96) =
    ## Parse .m4 file.
    ## Write only the overlaps which pass the filters.
    ## (no overlaps involving contained reads;
    ##  no overlaps involving short reads;
    ##  no overlaps with low identity)

    # Also, write a single letter for 5/3/C/O instead of "contained/contains/overlap".

    let sin = streams.newFileStream(in_fn, fmRead)
    defer: sin.close()
    let sout = streams.newFileStream(out_fn, fmWrite)
    defer: sout.close()

    let n = m4filtContainedStreams(sin, sout, min_len, min_idt_pct)

    # Write the number of lines into a separate file (to aid parsing, someday).
    let
        nout_fn = out_fn & ".n"
        fh = open(nout_fn, fmWrite)
    fh.write($n)
    fh.close()

proc split(index: M4Index, nProc: int): seq[int] =
    # Split index into at most nProc subsets, weighted by the
    # overlap count of each pile.
    # len(result) will be <= nProc

    var sizes: seq[int]
    for rec in index:
        sizes.add(rec.count)
    return util.splitWeighted(nProc, sizes)

proc m4filt(inFn: string,
 nProc: int,
 opts: M4filtOptions) =
    let index = getM4Index(inFn)
    let threadpool = threadpool_simple.newThreadPool(nProc)
    let splits = split(index, nProc)
    m4filtSingleton(inFn, splits, index, opts, threadpool)

proc m4filtRunner*(
 idtStage1: float = 90.0,
 idtStage2: float = 90.0,
 minLen: int = 6000,
 minCov: int = 2,
 maxCov: int = 200,
 maxDiff: int = 100,
 bestN: int = 10,
 minOverhang: int = 0,
 minDepth: int = 2,
 gapFilt: bool = false,
 keepIntermediates: bool = false,
 nProc: int = 24,
 inFn: string,
 filterLogFn: string,
 outFn: string) =
    ## Run the multi-stage m4 overlap filter (for HiFi Asm).
    ## Take only one m4 file, inFn. If in+'.idx' exists, that is the index.
    ## In stage one, reads that
    ## trigger a filter are marked including containment, gaps in coverage along the
    ## A-read, and repeat reads.
    ## In stage two the filters are applied and the N-best
    ## overlaps are kept for the 5prime and 3prime of each read.

    let opts = M4filtOptions(
        idtStage1: idtStage1,
        idtStage2: idtStage2,
        minLen: minLen,
        minCov: minCov,
        maxCov: maxCov,
        maxDiff: maxDiff,
        bestN: bestN,
        minOverhang: minOverhang,
        minDepth: minDepth,
        gapFilt: gapFilt,
        keepIntermediates: keepIntermediates,
        filterLogFn: filterLogFn,
        outFn: outFn)

    m4filt(inFn, nProc, opts)

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
 minOverhang: int = 0,
 gapFilt: bool = false,
 keepIntermediates: bool = false,
 nProc: int = 24,
 filterLogFn: string,
 debugLogFn: string = "LA4Falcon.log",
 outFn: string) =
    ##Runs the multi-stage m4 overlap filter for falcon. In stage one, reads that
    ##trigger a filter are marked including containment, gaps in coverage along the
    ##A-read, and repeat reads. In stage two the filters are applied and the N-best
    ##overlaps are kept for the 5prime and 3prime of each read.

    util.which("LA4Falcon")

    #LasJson is generated by the falcon pipeline
    let js = parseFile(lasJsonFn)
    var icmds: seq[string]
    for las in js:
        # Remember to redirect stderr, since we will not read it.
        # Otherwise, beyond several thousand lines, the stderr buffer
        # will fill up and block writing to stdout.
        let icmd = "LA4Falcon -mo {db} {las.getStr()} 2> {debugLogFn}".fmt
        icmds.add(icmd)
    let opts = M4filtOptions(
        idtStage1: idtStage1,
        idtStage2: idtStage2,
        minLen: minLen,
        minCov: minCov,
        maxCov: maxCov,
        maxDiff: maxDiff,
        bestN: bestN,
        minOverhang: minOverhang,
        minDepth: minDepth,
        gapFilt: gapFilt,
        keepIntermediates: keepIntermediates,
        filterLogFn: filterLogFn,
        outFn: outFn)

    let threadpool = threadpool_simple.newThreadPool(nProc)

    m4filtPiped(icmds, opts, threadpool)

proc ipaRunner*(ovlsFofnFn: string,
 idtStage1: float = 90.0,
 idtStage2: float = 90.0,
 minLen: int = 6000,
 minCov: int = 2,
 maxCov: int = 200,
 maxDiff: int = 100,
 bestN: int = 10,
 minOverhang: int = 0,
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
        minOverhang: minOverhang,
        minDepth: minDepth,
        gapFilt: gapFilt,
        keepIntermediates: keepIntermediates,
        filterLogFn: filterLogFn,
        outFn: outputFn)

    let threadpool = threadpool_simple.newThreadPool(nProc)

    m4filtPiped(icmds, opts, threadpool)

proc idx*(in_fn: string) =
    ## Given foo.m4, create index file foo.m4.idx
    ## (Or just use it if it exists.) Format is
    ## "0000 count start len", where count is the number of overlaps in the pile.
    ## (First field is dummy, instead of Aread name, for efficiency.)
    ## Return the sum of all pileups, which should exactly match filesize.
    let m4idx = op.getM4Index(in_fn)
    log("Number of piles={len(m4idx)}".fmt)

proc main() =
    echo "main"

when isMainModule:
    main()
