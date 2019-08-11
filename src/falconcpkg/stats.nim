# vim: sw=4 ts=4 sts=4 tw=0 et:
from ./util import nil
from algorithm import nil
from hts import `[]`
from json import `%*`
from math import nil
from sequtils import nil
from tables import `[]=`

proc get_all_ctgs(fasta: hts.Fai): tables.Table[string, int32] =
    ## Given indexed FASTA, return all name->length
    # (Hard to test, but very fast, since htslib can fseek to each header.)
    let n = hts.len(fasta)
    result = tables.initTable[string, int32]()
    for i in 0 ..< n:
        let name = fasta[i]
        result[name] = hts.chrom_len(fasta, name).int32

proc get_all_ctgs*(fasta_fn: string): tables.Table[string, int32] =
    var idx: hts.Fai
    if not hts.open(idx, fasta_fn):
        util.raiseEx("Could not open FASTA '" & fasta_fn & "'")
    return get_all_ctgs(idx)
# get_all_ctgs() were both copied from pbreports.

proc percentile(reverse_sorted_lens: seq[int32], p: float): int32 =
    let n = len(reverse_sorted_lens)
    var idx = (n.float64 * (1.0 - p)).int32
    if idx < 0: idx = 0
    if idx >= n: idx = n.int32 - 1
    return reverse_sorted_lens[idx]

proc esize(lens: seq[int32], total: var int64): float32 =
    var total: int64 = 0
    var sum_squares: int64 = 0
    for l in lens:
        let len = l.int64
        total += len
        sum_squares += len * len
    return sum_squares.float32 / total.float32

type
    NLStat* = tuple[nx: int32, lx: int]

proc len_stat*(reverse_sorted_lens: seq[int32], threshold: int32): NLStat =
    ## len_stat(reads, total * .50) would return (n50, l50)
    ## Remember: l50 is a count, not an index.
    var subtotal: int32 = 0
    for i in 0 ..< reverse_sorted_lens.len():
        let rl = reverse_sorted_lens[i]
        subtotal += rl
        if subtotal >= threshold:
            return (nx: rl, lx: i+1)
    util.raiseEx("threshold too high:" & $threshold & "; total=" & $subtotal)

type
    Stats = object
        # For now we ignore "N" bases, so we do not count gaps.
        max*, min*, median*, number*: int32
        sum*: int64
        mean*: float32
        nl50*: NLStat
        nl60: NLStat
        nl70: NLStat
        nl80: NLStat
        nl90: NLStat
        nl100*: NLStat
        esize: float32


proc to_json*(st: Stats): string =
    var j = %*
        {
            "sum": st.sum,
            "mean": math.round(st.mean),
            "median": st.median,
            "max": st.max,
            "N100": st.nl100.nx,
            "L100": st.nl100.lx,
            "N90": st.nl90.nx,
            "L90": st.nl90.lx,
            "N50": st.nl50.nx,
            "L50": st.nl50.lx,
            "esize": math.round(st.esize),
        }
    return json.pretty(j, indent=4)

proc calc_stats*(unsorted_lengths: seq[int32]): Stats =
    var lengths = algorithm.sorted(unsorted_lengths, order=algorithm.Descending)
    if lengths.len() == 0:
        return
    result.number = lengths.len().int32
    result.max = lengths[lengths.low.int32]
    result.min = lengths[lengths.high.int32]
    result.median = percentile(lengths, 0.50)
    result.sum = 0

    for length in lengths:
        result.sum += length
    result.mean = 1.0 * result.sum.float32 / result.number.float32;
    result.esize = esize(lengths, result.sum)

    var k = 1
    var cumulative: float64 = 0.0
    let sumf = result.sum.float64

    result.nl50 = len_stat(lengths, math.ceil(result.sum.float * 0.50).int32)
    result.nl60 = len_stat(lengths, math.ceil(result.sum.float * 0.60).int32)
    result.nl70 = len_stat(lengths, math.ceil(result.sum.float * 0.70).int32)
    result.nl80 = len_stat(lengths, math.ceil(result.sum.float * 0.80).int32)
    result.nl90 = len_stat(lengths, math.ceil(result.sum.float * 0.90).int32)
    result.nl100 = len_stat(lengths, math.ceil(result.sum.float * 1.00).int32)

proc assembly*(fasta_fn: string) =
    let name2length = get_all_ctgs(fasta_fn)
    var lengths = sequtils.toSeq(tables.values(name2length))
    let stats = calc_stats(lengths)
    echo to_json(stats)
