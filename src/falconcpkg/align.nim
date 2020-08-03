# vim: sw=4 ts=4 sts=4 tw=0 et:
import hts
from hts/private/hts_concat import bam_cigar2qlen, bam_cigar2rlen
from algorithm import nil
from parseutils import nil
from sequtils import nil
from strutils import format
from tables import contains, `[]`, `[]=`
from sets import items
from streams import nil
from strformat import fmt
from ./util import log, raiseEx

proc logRec(record: Record) =
    # I think len == stop-start+1, but I need to verify. ~cd
    var s: string
    discard hts.sequence(record, s)
    log(format("$# $# ($#) [$# .. $#] seqlen=$# $#...", record.tid,
            record.chrom, record.qname, record.start, record.stop, s.len(),
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


proc target_to_query*(cigar: hts.Cigar, t_pos, t_start, t_end: int64): tuple[
        q_start, q_end: int64] =
    #given a region in target space return the coordinates of the query for the same region.
    #for more info see page seven of the bible :  https://samtools.github.io/hts-specs/SAMv1.pdf

    var t_off, q_off: int64 = 0
    t_off += t_pos
    result.q_start = 0
    var start_unset, end_unset: bool = true
    var tdiffs, tdiffe: int64 = 0
    var ncig: int = 0
    var last: CigarElement
    for op in cigar:
        last = op
        inc(ncig)
        let cons = op.consumes
        if cons.query:
            q_off += op.len
        if cons.reference:
            t_off += op.len
        tdiffs = t_off - t_start
        tdiffe = t_end - t_off
        if t_off >= t_start and start_unset:
            if cons.query and cons.reference:
                if t_pos < t_start:
                    #echo "A t_start:", t_start, "t_off:", t_off, " tdiffs:",  tdiffs , " qoff:", q_off, " op:", op, " ncig:", ncig
                    result.q_start = qoff - tdiffs
            else:
                result.q_start = q_off
            start_unset = false
        if t_off >= t_end and end_unset:
            if cons.query and cons.reference:
                #echo "X t_start:", t_start, "t_off:", t_off, " tdiffe:",  tdiffe , " qoff:", q_off, " op:", op, " ncig:", ncig
                result.q_end = (q_off + tdiffe) - 1
            else:
                result.q_end = q_off

            end_unset = false
            break

    if end_unset:
        if last.consumes.reference:
            result.q_end = q_off - 1
        else:
            result.q_end = (q_off - last.len) - 1

type
    # based on https://github.com/zeeev/bamPals/blob/master/src/enrich_optional_tags.c
    Pal = object
        ref_beg, ref_end, ref_len: int32 # alignment in reference coordinates
        qry_len: int32                   # "I" tag in BAM
        pct_idt: float32                 # "f" tag in BAM

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
    result.ref_beg = hts.start(record).int32
    result.ref_end = hts.stop(record).int32
    result.ref_len = result.ref_end - result.ref_beg
    result.qry_len = bam_cigar2qlen(record.b.core.n_cigar.cint,
            hts.bam_get_cigar(record.b))

    let xr = bam_cigar2rlen(record.b.core.n_cigar.cint, hts.bam_get_cigar(
            record.b))
    assert result.ref_len == xr

    let (t_consumed, q_sclipped) = pal_cigar(record.cigar)
    assert result.ref_len == t_consumed

    let nmtag = hts.tag[int](record, "NM")
    if hts.isSome(nmtag):
        let edit_dist = hts.get(nmtag)
        let pi: float = 100.0 *
            (record.b.core.l_qseq - q_sclipped - edit_dist).float /
            (record.b.core.l_qseq - q_sclipped).float
        result.pct_idt = pi
    else:
        result.pct_idt = -1.0

proc tags_enrich(record: hts.Record) =
    # Mutate the underlying object.

    #let (qstart, qend, qlen) = calc_query_pos(record)
    #log(format(" calc(qstart=$#, qend=$#, qlen=$#)", qstart, qend, qlen))
    let pal = pal_calc(record)
    hts.set_tag[int](record, "XB", pal.ref_beg) # same as core.pos
    hts.set_tag[int](record, "XE", pal.ref_end) # same as bam_endpos
    hts.set_tag[int](record, "XR", pal.ref_len) # same as bam_cigar2rlen
    hts.set_tag[int](record, "XQ", pal.qry_len) # same as bam_cigar2qlen
    if pal.pct_idt >= 0:
        hts.set_tag[float](record, "XP", pal.pct_idt)

proc bam_tags_enrich(new_bam: var hts.Bam, old_bam: hts.Bam) =
    hts.write_header(new_bam, old_bam.hdr)

    for record in old_bam:
        tags_enrich(record)
        hts.write(new_bam, record)

proc bam_tags_enrich*(output_fn, input_fn: string) =
    ## Add XB/XE/XP/XR/XQ: beg/end/%idt/aln-ref-len/qry-len
    var
        obam, ibam: hts.Bam
    hts.open(obam, output_fn, mode = "w") # compression will be the default for the format
    hts.open(ibam, input_fn)

    try:
        bam_tags_enrich(obam, ibam)
    finally:
        hts.close(ibam)
        hts.close(obam)

let
    SOFT_CLIP = 'S'
    HARD_CLIP = 'H'
    CLIPPINGOPS = [SOFT_CLIP, HARD_CLIP]

proc get_left_right_clip(cigar: hts.Cigar): tuple[left, right: int] =
    var n = hts.len(cigar)

    var leftclip = 0
    for i in countup(0, n - 1):
        let
            elem = cigar[i]
            op = $op(elem)
            opLen = len(elem)
        if op notin CLIPPINGOPS:
            break
        leftclip += oplen

    var rightclip = 0
    for i in countdown(n - 1, 0):
        let
            elem = cigar[i]
            op = $op(elem)
            opLen = len(elem)
        if op notin CLIPPINGOPS:
            break
        rightclip += oplen

    return (left: leftclip, right: rightclip)

proc get_reference_length(record: hts.Record, targets: seq[Target]): uint32 =
    let
        tid = hts.tid(record)
        reference_target = targets[tid]
    return reference_target.length

proc bam_count*(input_fn: string): int =
    ## Simply return the number of records.
    var
        in_bam: hts.Bam
        n: int
    hts.open(in_bam, input_fn)
    defer: hts.close(in_bam)

    for record in in_bam:
        n += 1
    return n

proc bam_filter_clipped(ocount: File, obam: var hts.Bam, ibam: hts.Bam,
        max_clipping, end_margin: int, flags_exclude: uint16, verbose,
                tags_enrich: bool): int =
    # Return the number skipped. Write the rest into obam.
    hts.write_header(obam, ibam.hdr)

    let targets = hts.targets(ibam.hdr)

    if verbose:
        log(format(
                "Listing skipped records, with max_clipping=$# and end_margin=$# ...",
            max_clipping, end_margin))
        log(
                "tid chrom (qname) [start .. end+1 (0-based)] seqlen cigar(truncated)")
    var n_skipped = 0
    var n_kept = 0

    for record in ibam:
        let flag = hts.flag(record)
        if hts.has_flag(flag, flags_exclude):
            if verbose:
                logRec(record)
            n_skipped += 1
            continue

        let (leftclip, rightclip) = get_left_right_clip(record.cigar)

        let pal = pal_calc(record)

        if leftclip > max_clipping:
            # left clipping is high
            if pal.ref_beg > end_margin:
                # alignment left is not near the contig start
                if verbose:
                    logRec(record)
                n_skipped += 1
                continue

        if rightclip > max_clipping:
            # right clipping is high
            let reference_length = get_reference_length(record, targets).int
            if (pal.ref_end + end_margin) < reference_length:
                # alignment right is not near the contig end
                if verbose:
                    logRec(record)
                n_skipped += 1
                continue

        if tags_enrich:
            tags_enrich(record)
        inc(n_kept)
        hts.write(obam, record)

    ocount.writeline(n_kept)
    return n_skipped

proc toUInt16*(v: string): uint16 =
    var number: int
    if strutils.startsWith(v, "0x"):
        let digits = parseutils.parseHex(v, number, 2)
    else:
        let digits = parseutils.parseInt(v, number)
    return number.uint16

proc bam_filter_clipped*(output_count_fn, output_fn, input_fn: string,
        max_clipping = 100, end_margin = 25, Flags_exclude = "0",
                verbose = false, tags_enrich = false) =
    ## Filter alignments with significant clipping.
    ## To skip an alignment, both max_clipping and end_margin must be exceeded on at least 1 end.
    let
        flags_exclude: uint16 = toUInt16(Flags_exclude)
    var
        obam, ibam: hts.Bam
    hts.open(obam, output_fn, mode = "w") # compression will be the default for the format
    hts.open(ibam, input_fn)
    let ocount = open(output_count_fn, fmWrite)

    try:
        let n_skipped = bam_filter_clipped(ocount, obam, ibam,
            max_clipping, end_margin, flags_exclude, verbose, tags_enrich)
        if verbose:
            log("Skipped:", n_skipped)
    finally:
        hts.close(ibam)
        hts.close(obam)
        ocount.close()


# Someday we might want to look at the CIGAR alignments.
# But for now, we call it a decent alignment if it is long enough.
proc update_counts(bam_fn: string, params: Params,
        mapped: var sets.HashSet[string],
        exclusions: var tables.CountTable[string]) =
    var
        totals: tables.Table[string, int]
        qlengths: tables.Table[string, int]
        in_bam: hts.Bam
    totals = tables.initTable[string, int]()
    qlengths = tables.initTable[string, int]()

    hts.open(in_bam, bam_fn)
    defer: hts.close(in_bam)

    for record in in_bam:
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

proc bam_filter_ipa*(bams_fofn: string, all_subread_names_fn = "",
        min_len = -1,
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

type
    PafLine = object
        qname: string
        qlen, qstart, qend: int
        relStrand: char
        tname: string
        tlen, tstart, tend: int
        numResidueMatches: int
        alnBlockLen: int
        mappingQuality: int

# https://bioconvert.readthedocs.io/en/master/formats.html#paf-pairwise-mapping-format
proc `$`(p: PafLine): string =
    return "{p.qname}\t{p.qlen}\t{p.qstart}\t{p.qend}\t{p.relStrand}\t{p.tname}\t{p.tlen}\t{p.tstart}\t{p.tend}\t{p.numResidueMatches}\t{p.alnBlockLen}\t{p.mappingQuality}".fmt

proc writePaf(out_paf: streams.Stream, record: hts.Record, targets: seq[Target]) =
    var p: PafLine

    let qpos = calc_query_pos(record)
    p.qname = record.qname
    p.qstart = qpos.qstart
    p.qend = qpos.qend
    p.qlen = qpos.qlen # same as bam_cigar2qlen()
    p.relStrand = '+'

    if record.flag.reverse: # reversed
        let clipleft = p.qstart
        let clipright = p.qlen - p.qend
        p.qstart = clipright
        p.qend = p.qlen - p.qstart
        p.relStrand = '-'

    p.tname = record.chrom
    p.tstart = hts.start(record).int
    p.tend = hts.stop(record).int
    #p.tlen = p.tend - p.tstart  # not what PAF wants
    let reference_length = get_reference_length(record, targets).int
    p.tlen = reference_length

    p.alnBlockLen = p.tend - p.tstart # same as bam_cigar2rlen()

    var mismatches = 0
    if not hts.isNone(hts.tag[int](record, "NM")):
        mismatches = hts.tag[int](record, "NM").get.int
    else:
        mismatches = 0 # not correct, but also not important
    p.numResidueMatches = p.alnBlockLen - mismatches
    p.mappingQuality = record.mapping_quality.int
    streams.write(out_paf, $p)
    streams.write(out_paf, '\n')

proc bam2paf*(in_bam_fn, out_p_paf_fn: string, out_a_paf_fn: string) =
    ## https://bioconvert.readthedocs.io/en/master/formats.html#paf-pairwise-mapping-format
    var b: Bam
    open(b, in_bam_fn, index = true)
    defer:
        b.close()
    var out_p_paf = streams.openFileStream(out_p_paf_fn, fmWrite)
    var out_a_paf = streams.openFileStream(out_a_paf_fn, fmWrite)
    defer:
        streams.close(out_p_paf)
        streams.close(out_a_paf)
    let targets = hts.targets(b.hdr) # in case there are more than 1
    for record in b:
        let rname = record.chrom
        if strutils.find(rname, '-') >= 0:
            writePaf(out_a_paf, record, targets)
        else:
            writePaf(out_p_paf, record, targets)
    # This is fine, but canonical PAF would add these tags also:
    # var extra = ["mm:i:"+(NM-I[1]-D[1]), "io:i:"+I[0], "in:i:"+I[1], "do:i:"+D[0], "dn:i:"+D[1]];
    # https://github.com/lh3/miniasm/blob/master/misc/sam2paf.js
