# vim: sw=4 ts=4 sts=4 tw=0 et:
import hts
from algorithm import nil
from sequtils import nil
from strutils import format
from tables import contains, `[]`, `[]=`
from sets import items
from ./util import log

proc logRec(record: Record) =
    # I think len == stop-start+1, but I need to verify. ~cd
    var s: string
    discard hts.sequence(record, s)
    log(format("$# $# ($#) [$# .. $#] $# seqlen=$#", record.tid, record.chrom,
            record.qname, record.start, record.stop,
        ($record.cigar).substr(0, 32), s.len()))

type
    Params = object
        min_frac: float

## At each end of Cigar, we return
##  * +N for N soft-clip
##  * -N for N hard-clip
##  *  0 for other Cigar (or Cigar.len < 2)
proc clips*(cigar: hts.Cigar): tuple[left: int, right: int] =
    var left, right: int
    left = 0
    right = 0

    if (0 == cigar.len): return (left, right)

    let first: hts.CigarElement = cigar[0]
    if 'H' == $first.op:
        left = -first.len
    elif 'S' == $first.op:
        left = first.len

    if (1 == cigar.len): return (left, right)

    let second: hts.CigarElement = cigar[cigar.len - 1]
    if 'H' == $second.op:
        right = -second.len
    elif 'S' == $second.op:
        right = second.len

    return (left, right)

proc calc_query_pos*(record: hts.Record): tuple[qstart: int, qend: int,
        qlen: int] =
    var q: string # temp
    discard hts.sequence(record, q)

    var qlen: int = len(q)
    var qstart: int = 0
    var qend: int = qlen
    let (q_clip_left, q_clip_right) = clips(record.cigar)

    if q_clip_left < 0:
        qstart -= q_clip_left
        qlen -= q_clip_left
        qend -= q_clip_left
    else:
        qstart += q_clip_left

    if q_clip_right < 0:
        qlen -= q_clip_right
    else:
        qend -= q_clip_right

    # echo "q_clip_left = ", q_clip_left, ", q_clip_right = ", q_clip_right, ", len(q) = ", len(q)
    # echo "qstart = ", qstart, ", qend = ", qend, ", qlen = ", qlen

    return (qstart, qend, qlen)

# Someday we might want to look at the CIGAR alignments.
# But for now, we call it a decent alignment if it is long enough.
proc update_counts(bam_fn: string, params: Params,
        mapped: var sets.HashSet[string],
        exclusions: var tables.CountTable[string]) =
    var
        totals: tables.Table[string, int]
        qlengths: tables.Table[string, int]
        b: hts.Bam
    totals = tables.initTable[string, int]()
    qlengths = tables.initTable[string, int]()

    hts.open(b, bam_fn)
    defer: hts.close(b)
    for record in b:
        let key: string = record.qname
        let (qstart, qend, qlen) = calc_query_pos(record)

        sets.incl(mapped, key)

        # let refspan = record.stop - record.start + 1
        let qspan = qend - qstart
        totals[key] = qspan + tables.getOrDefault(totals, key, 0)
        qlengths[key] = qlen
        #logRec(record)
    for key, qlen in tables.pairs(qlengths):
        let span_total = totals[key]
        let frac = float(span_total) / float(qlen)
        if frac >= params.min_frac:
            tables.inc(exclusions, key)
            #log("exclude:" & key)
        #log(format("$# $#/$#=$#", key, reflentotal, qlen, frac))

    #proc find_all(fn: string): sets.HashSet[string] =
    #    #result = sets.initHashSet[string]()
    #    if fn == "":
    #        return
    #    log(format("Reading from '$#'", fn))
    #    for line in lines(fn):
    #        sets.incl(result, line)

proc align_filter*(bams_fofn: string, all_subread_names_fn = "", min_len = -1,
        min_frac = 0.70) =
    ## Print subreads which have decent alignments in any of the bam inputs.
    if min_len != -1:
        log("WARNING: --min-len ignored")
    var params = Params(min_frac: min_frac)
    log(params)
    #var all_subread_names = find_all(all_subread_names_fn)
    #log("all_subread_names.len=", sets.len(all_subread_names))
    var mapped = sets.initHashSet[string]()
    var exclusions = tables.initCountTable[string]()
    for fn in lines(bams_fofn):
        log("Processing ", fn)
        update_counts(fn, params, mapped, exclusions)
    log("mapped.len=", sets.len(mapped))
    log("exclusions.len=", tables.len(exclusions))
    #var unmapped = sets.difference(all_subread_names, mapped)
    #log("unmapped.len=", sets.len(mapped))
    #for key in sets.items(unmapped):
    #    tables.inc(exclusions, key)
    #log("exclusions.len=", tables.len(exclusions))
    var sorted_exclusions: seq[string] = sequtils.toSeq(tables.keys(
            exclusions))
    algorithm.sort(sorted_exclusions)
    for key in sorted_exclusions:
        #log("Count of '", key, "': ", everything[key])
        echo key
