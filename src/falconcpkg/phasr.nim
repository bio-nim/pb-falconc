# vim: sw=4 ts=4 sts=4 tw=0 et:
from algorithm import nil
from sequtils import nil
from strutils import format
from strformat import fmt
import sets
import tables
import deques
import hts
#import os
import random
import math

from ./util import raiseEx
import networkx/classes/wgraph
import networkx/classes/graph
from networkx/algorithms/components/connectedc import connected
from "./kmers" import pot_t, spot_t, Dna, dna_to_kmers, initSpot, difference,
    uniqueShared, nuniq #import kmers

type
    opts = object
        kmersize, iterations: int
        delta: float
        out_fn: string

var globalOpts = opts()

var inverse = newTable[int, string]() # readid -> name
var rpos = newTable[string, int]()

proc getWithin(u: Node, g: ref Graph[int], phase: TableRef[int, int]): float =
    var
        within, totalc, totalw: int = 0
        phases = initCountTable[int]()

    for v, w in successors(g, u):
        inc(totalc)
        phases.inc(phase[v])
        if phase[u] == phase[v]:
            within += w
        totalw += w

    result = float(within)/float(totalw)


proc ranScore(n: int, scores: var seq[float]): float =
    result = 0
    shuffle(scores)
    for i in 1..n:
        result += scores[i]
    result = result / float(n)

proc switchState(n: Node, g: ref Graph[int], p: TableRef[int, int],
        ploidy: int, iter: int, temp: float) =

    var
        phaseW = initTable[int, int]()
        phaseC = initTable[int, int]()
        phaseD = initTable[int, float]()
        denSum: float = 0

    for i in 0..ploidy-1:
        phaseW[i] = 0
        phaseC[i] = 0
        phaseD[i] = 0.0

    for v, w in successors(g, n):
        phaseW[p[v]] += w
        phaseC[p[v]] += 1

    for i in 0..ploidy-1:
        phaseD[i] = float(phaseW[i]) / float(phaseC[i])
        denSum += phaseD[i]

    for i in 0..ploidy-1:
        phaseD[i] /= denSum

    var newPhase = (p[n] + 1 + rand(ploidy - 2)) mod ploidy
    var prop: float = pow(E, ((phaseD[newPhase] - phaseD[p[n]])/temp))

    if iter mod 1000 == 0:
        echo prop, " ", phaseD[newPhase], " ", " ", phaseD[p[n]], " ", temp,
            " ", iter

    if prop > rand(1.0):
        p[n] = newPhase



proc overallScore(g: ref Graph[int], p: TableRef[int, int]): float =
    var overall: float = 0
    for n in nodes(g):
        #echo "n:", n
        overall += getWithin(n, g, p)
        #echo "overall:", overall
    return overall

type
    PhaseStats = object
        hist: seq[int] # time spent in each phase
        total: int     # total time

proc algo(
    rn: TableRef[string, int], # rn:  name -> readid
    g: ref Graph[int] # Node is readid
    ): Table[int, int] = # phase-block-id -> local-phase-id (up to ploidy)


    var phase = newTable[int, int]() # Node -> local-phase-id
    var stats = newTable[int, PhaseStats]()             # Node -> local-phase-id

    const ploidy = 2

    for readname, i in rn:
        phase[i] = rand(ploidy-1) # up to ploidy
        stats[i] = PhaseStats(hist: newSeq[int](ploidy))

    var nSeq = newSeq[Node]()
    for v in nodes(g):
        nSeq.add(v)

    var
        temp, minTemp, alpha: float = 0
    temp = 1.0
    alpha = 0.90
    minTemp = 0.0001


    for i in 1..globalOpts.iterations:
        shuffle(nSeq)
        for v in nSeq:
            switchState(v, g, phase, ploidy = ploidy, i, temp)
            stats[v].total += 1
            stats[v].hist[phase[v]] += 1
        #  echo "PD\t{i}\t{v}\t{phase[v]}\t{inverse[v]}".fmt
        var score = overallScore(g, phase)
        if (i mod 100) == 0:
            echo fmt("currently at i={i}, score={score:0.3f}")
            temp *= alpha
            temp = max(minTemp, temp)

    var output = open(globalOpts.out_fn, fmWrite)
    result = initTable[int, int]()
    var frqRes = initTable[int, float]()
    var histf = newSeq[float](2) # assume ploidy==2 for now
    for rid in algorithm.sorted(sequtils.toSeq(values(rn))):
        var max: float = -1.0
        for i in 0 ..< ploidy:
            histf[i] = stats[rid].hist[i] / stats[rid].total
            if histf[i] > max:
                result[rid] = i
                max = histf[i]
                frqRes[rid] = histf[i]
        # Read was removed before algorithm (zero weights across all edges). This shoudl not be hit, but will remain for satefy.
        if not (rid in frqRes):
            output.write("-1\t-1\t-1\t\t0.0\t{inverse[rid]}\n".fmt)
            continue
        # If the Frequency is not high enough node is removed.
        if frqRes[rid] < globalOpts.delta:
            output.write("-1\t-1\t-1\t\t{frqRes[rid]}\t{inverse[rid]}\n".fmt)
            wgraph.remove_node(g, rid.Node)

    var seen = sets.initHashSet[int]()
    var phaseBlock = 0
    for v in nodes(g):
        if seen.contains(v):
            continue
        inc(phaseBlock)
        for n in connected(g, v):
            if seen.contains(n):
                continue
            seen.incl(n)
            output.write(
                    "{n}\t{phaseBlock}\t{result[n]}\t{frqRes[n]}\t{inverse[n]}\n".fmt)
    output.close()

proc showRec(record: Record) =
    echo format("$# $# ($#) [$# .. $#] $#", record.tid, record.chrom,
            record.qname, record.start, record.stop,
        ($record.cigar).substr(0, 32))
type
    ProcessedRecord = object
        # These are both refs already.
        rec: Record
        kmers: kmers.spot_t

proc processRecord(record: Record, klen: int, rseq: Dna): ProcessedRecord =
    var qseq: string
    discard hts.sequence(record, qseq)
    let kmers: pot_t = dna_to_kmers(qseq, klen)
    #echo "kmers.seeds.len=", kmers.seeds.len(), " qseq.len=", qseq.len(), " aln.len", (record.stop-record.start)
    # Note that record.stop is 1 past the last index.

    # complement to remove kmers that are in the reference where the read maps.
    var rsubseq = rseq.substr(record.start - klen + 1,
            record.stop-1 + klen - 1) # TODO(CD): Use a slice?
                                      #var rsubseq = rseq.substr(record.start, record.stop-1)  # TODO(CD): Use a slice?
    var refkmers = dna_to_kmers(rsubseq, klen)
    #echo "refkmers.len=", refkmers.seeds.len(), "(", rsubseq.len(), "), kmers.len=", kmers.seeds.len(), "(", qseq.len(), ")", (record.stop-record.start)
    var refspot = initSpot(refkmers)
    var filtKmers = difference(kmers, refspot)
    #[echo " Filtered kmers.len=", filtKmers.seeds.len(), " from ",
        kmers.seeds.len(), " for read ", record.qname
        ]#
    var finalKmers = initSpot(filtKmers)
    return ProcessedRecord(rec: record, kmers: finalKmers)

type
    Pileup = object
        # I'm not sure I've named this well. It's a list of records that overlap a given one,
        # all processed for kmers already.
        precs: seq[ProcessedRecord]
        current: ProcessedRecord
        # All other recs overlap this one.
        min_start, max_stop: int

proc filterPileup(queue: Deque[ProcessedRecord], cqi: int): Pileup =
    # We know nothing is completely to the left of current because all those
    # are popped after each "yield". However, recs could be to right of current,
    # although they overlapped something longer, earlier in the queue.
    let current = queue[cqi]
    let current_stop = current.rec.stop
    result.precs = newSeqOfCap[ProcessedRecord](queue.len())
    result.current = current
    var
        min_start = current.rec.stop
        max_stop = current.rec.start
    for pr in queue.items():
        if pr == current or pr.rec.start >= current_stop:
            continue # must overlap on the right
        if min_start > pr.rec.start:
            min_start = pr.rec.start
        if max_stop < pr.rec.stop:
            max_stop = pr.rec.stop
        result.precs.add(pr)
    result.max_stop = max_stop
    result.min_start = min_start
    assert result.min_start == result.precs[0].rec.start # since we popped the rest already

type
    RecordGetter = iterator (b: hts.Bam): hts.Record

iterator nextRecord(b: hts.Bam): hts.Record {.closure.} =
    for r in hts.items(b): # Note: hts.items is not a closure iterator.
        #TODO allow user to specify
        if (r.mapping_quality < 20):
            continue
        if (r.flag.int and 3844) > 0:
            continue
        if r.start >= r.stop:
            continue # or raise exception?
                 # Without this check, zero-length records could cause problems later.
                 # Note: In hts, "stop" is one past the last position.
        rpos[r.qname] = r.start
        yield r

proc nextBamRecord(next: RecordGetter, b: hts.Bam, r: var hts.Record): bool =
    if system.finished(next):
        return false
    let nr = next(b)
    if system.finished(next):
        return false
    r = hts.copy(nr) # need a copy because iterator stores record on the Bam struct
    #showRec(r)  # DEBUG
    #echo "new_rec:", r.qname
    return true

iterator overlaps(b: hts.Bam, klen: int, rseq: string): Pileup =
    var next: RecordGetter = nextRecord # Instantiate closure iterator.
    var last_stop = 0
    var queue = deques.initDeque[ProcessedRecord](64)
    var current_queue_index: int = 0 # initially not real

    while true:
        #  Ensure that current_queue_index points at a real record,
        #  or quit.
        if queue.len() <= current_queue_index:
            var new_record: hts.Record = nil
            let gotAnother = nextBamRecord(next, b, new_record)
            if not gotAnother:
                break # None left!
            queue.addLast(processRecord(new_record, globalOpts.kmersize, rseq))
        assert queue.len() > current_queue_index
        let current = queue[current_queue_index].rec

        # Flush (pop from left) any records which do not overlap current.
        while queue.peekFirst().rec.stop <=
                current.start: # and assume current is not zero-length
            discard queue.popFirst()
            current_queue_index -= 1
            assert current == queue[current_queue_index].rec

        #echo " first: queue.len=", queue.len()
        #showRec(queue.peekFirst().rec)
        #echo " curr: current_queue_index=", current_queue_index
        #showRec(current)
        #assert current == queue[current_queue_index].rec
        #echo " last:"
        #showRec(queue.peekLast().rec)

        # Add records to end-of-queue until:
        # - last record starts later than current record stops, or
        # - we have no more records from Bam.
        while true:
            if queue.len() == 0 or queue.peekFirst().rec.start < current.stop:
                var new_record: hts.Record = nil
                let gotAnother = nextBamRecord(next, b, new_record)
                if not gotAnother:
                    break # but keep processing the queue
                queue.addLast(processRecord(new_record, globalOpts.kmersize,
                        rseq))

        # Yield the current record.
        #echo fmt("Yield at [{current.start}, {current.stop}) for current_queue_index={current_queue_index}")
        yield filterPileup(queue, current_queue_index) # YIELD
        current_queue_index += 1

## Count the number of kmers in both a and b.
## (Should we count only positional matches?)
#
proc jaccardDist*(a, b: ProcessedRecord): float =
    var num: int = uniqueShared(a.kmers, b.kmers)
    var da: int = nuniq(a.kmers)
    var db: int = nuniq(b.kmers)
    var den = da + db
    #echo format(" num/den=$#/($#+$#) == $#", num, da, db, num/den)
    return num/den

proc readaln*(bfn: string; fasta: string): tuple[t: TableRef[string, int],
        g: ref Graph[int]] =
    var b: hts.Bam

    let g = newGraph[int]()
    let readIdLookup = newTable[string, int]()
    var readCount: int = 1

    hts.open(b, bfn, index = true)
    # We do not really need the index, but it proves that the Bam is sorted.
    # Actually, we expect it to be sorted with @SO=coordinate, so maybe we
    # should verify that. TODO(CD).
    # Note that the sort is wrong for circular genomes; the secondary key is
    # the query, not the coordinate. So we will need to do something tricky
    # for circular references. TODO(CD).
    echo "[INFO] reading bam"
    for ovlps in overlaps(b, klen = globalOpts.kmersize, fasta.Dna):
        #echo fmt("ovlps.current.rec.qname='{ovlps.current.rec.qname}'")
        if not readIdLookup.hasKey(ovlps.current.rec.qname):
            readIdLookup[ovlps.current.rec.qname] = readCount
            inc(readCount)
        #echo format(" Processing pileup of $# precs", ovlps.precs.len())
        for i in 0 ..< ovlps.precs.len():
            #echo fmt(" ovlps.precs[{i}].rec.qname='{ovlps.precs[i].rec.qname}'")
            if not readIdLookup.hasKey(ovlps.precs[i].rec.qname):
                readIdLookup[ovlps.precs[i].rec.qname] = readCount
                inc(readCount)
            let jd = jaccardDist(ovlps.current, ovlps.precs[i])
            let rcount = int(jd * 1000) + 1

            echo fmt(
                    " Adding edge ( {ovlps.current.rec.qname} -> {ovlps.precs[i].rec.qname} ) w= {rcount}")

            wgraph.add_edge(g, rcount, (readIdLookup[
                    ovlps.current.rec.qname].Node,
                    readIdLookup[ovlps.precs[i].rec.qname].Node))
    for readname, i in readIdLookup:
        inverse[i] = readname
    result.t = readIdLookup
    result.g = g



proc main*(aln_fn: string, ref_fn: string, out_fn: string, iterations = 1000,
        kmersize = 21, delta = 0.75) =
    ##Phase PacBio CCS/HIFI reads.

    globalOpts.iterations = iterations
    globalOpts.kmersize = kmersize
    globalOpts.delta = delta
    globalOpts.out_fn = out_fn

    #prime random
    randomize()


    echo "[INFO] input reference (fasta):", ref_fn,
        "\n[INFO] alignment (reads):", aln_fn
    if strutils.find(ref_fn, "fa") == -1:
        echo format("[WARN] Bad fasta filename? '$#'", ref_fn)
    var refx: hts.Fai
    if not hts.open(refx, ref_fn):
        raiseEx(format("[FATAL] Could not open '$#'", ref_fn))
    #assert refx.len() == 1
    let reference_dna = refx.get(refx[0])
    var gdat = readaln(aln_fn, reference_dna)
    discard algo(gdat.t, gdat.g)
    echo "[INFO] all good, bye!"

when isMainModule:
    main()
