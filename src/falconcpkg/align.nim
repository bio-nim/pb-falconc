# vim: sw=4 ts=4 sts=4 tw=0 et:
import hts
from hts/private/hts_concat import bam_cigar2qlen, bam_cigar2rlen
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
    log(format("$# $# ($#) [$# .. $#] seqlen=$# $#...", record.tid, record.chrom,
            record.qname, record.start, record.stop, s.len(),
        ($record.cigar).substr(0, 32)))

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
    ## qlen is original query, so (qend-start) < qlen if clipped (hard or soft)
    ## seq excludes hard-clip, so len(seq) < qlen if hard-clipped
    ## (i.e. qlen is not well-named)
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

type
    # based on https://github.com/zeeev/bamPals/blob/master/src/enrich_optional_tags.c
    Pal = object
        xb, xe, xl: int32  # "I" tags in BAM
        xp: float32  # "f" tag in BAM

proc pal_cigar(cigar: hts.Cigar): tuple[t_consumed, q_sclipped: int32] =
    var t_consumed, q_sclipped: int
    for elem in cigar:
        case $hts.op(elem)
        of '=', 'D', 'X', 'M', 'N':
            t_consumed += hts.len(elem)
        of 'S':
            q_sclipped += hts.len(elem)
        else:
            discard
    return (t_consumed: t_consumed.int32, q_sclipped: q_sclipped.int32)

proc pal_calc(record: hts.Record): Pal =
    let (t_consumed, q_sclipped) = pal_cigar(record.cigar)
    result.xl = t_consumed
    result.xb = record.b.core.pos
    result.xe = result.xb + result.xl

    let nmtag = hts.tag[int](record, "NM")
    if hts.isSome(nmtag):
        let edit_dist = hts.get(nmtag)
        let pi: float = 100.0 *
            (record.b.core.l_qseq - q_sclipped - edit_dist).float /
            (record.b.core.l_qseq - q_sclipped).float
        result.xp = pi
    else:
        result.xp = -1.0
    
proc bam_tags_enrich*(old_bam_fn, new_bam_fn: string) =
    ## Add XB/XE/XP/XR: beg/end/%idt/aln-ref-len
    var
        ob, nb: hts.Bam
    hts.open(ob, old_bam_fn)
    hts.open(nb, new_bam_fn, mode="w")
    defer: hts.close(nb)
    defer: hts.close(ob)
    hts.write_header(nb, ob.hdr)

    for record in ob:
        # I want to remember how to use bam_cigar2qlen
        #logRec(record)
        #let (qstart, qend, qlen) = calc_query_pos(record)
        #log(format(" calc(qstart=$#, qend=$#, qlen=$#)", qstart, qend, qlen))
        #log(format(" cigar2qlen=$#", bam_cigar2qlen(record.b.core.n_cigar.cint, hts.bam_get_cigar(record.b))))
        #log(format(" cigar2rlen=$#", bam_cigar2rlen(record.b.core.n_cigar.cint, hts.bam_get_cigar(record.b))))
        let pal = pal_calc(record)
        #log(format(" pal XB=$# XE=$# XP=$# XL=$#", pal.xb, pal.xe, pal.xp, pal.xl))
        hts.set_tag[int](record, "XB", pal.xb) # same as core.pos
        hts.set_tag[int](record, "XE", pal.xe) # same as bam_endpos
        hts.set_tag[int](record, "XR", pal.xl) # same as bam_cigar2rlen?
        #hts.set_tag[int](record, "XQ", pal.xl) # same as bam_cigar2qlen
        if pal.xp >= 0:
            hts.set_tag[float](record, "XP", pal.xp)
        hts.write(nb, record)


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
