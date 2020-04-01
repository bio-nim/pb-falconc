# vim: sw=4 ts=4 sts=4 tw=0 et:
from algorithm import nil
from json import `%*`
from math import nil
import strutils
import tables
from stats import nil
import streams
import strformat

type
    Stats* = object
        max*, min*, median*, p95*, number*: int32
        sum*: int64
        mean*: float32
        nl50*: stats.NLStat
        nl95: stats.NLStat
        esize: float32
        coverage: float32

type
    SeqLenTuple* = tuple[id: int64, name: string, bases: int32]

proc parse_seq_len_tuples*(sin: Stream): seq[SeqLenTuple] =
    var ret = newSeq[SeqLenTuple]()
    var line = ""
    while sin.readLine(line):
        if line.len() == 0:
            continue
        if line[0] != 'S':
            continue
        let tokens = line.split('\t')
        if tokens.len() != 7:
            raise newException(Exception, fmt"Malformed sequence ('S') line in the RaptorDB stream. Line = '{line}'")
        ret.add( (id: parseInt(tokens[1]).int64, name: tokens[2], bases: parseInt(tokens[3]).int32))
    return ret

proc get_lens_from_seq_tuples*(seq_tuples: seq[SeqLenTuple]): seq[int32] =
    var ret = newSeq[int32]()
    for vals in seq_tuples:
        ret.add(vals[2]) # Add the bases.
    return ret

proc cutoff_lengths*(length_cutoff: int32, seq_lens: seq[int32]): seq[int32] =
    var ret = newSeq[int32]()
    for val in seq_lens:
        if val >= length_cutoff:
            ret.add(val)
    return ret

proc compute_truncation*(seq_tuples_raw, seq_tuples_preads: seq[SeqLenTuple]): float32 =
    # Sanity check.
    if seq_tuples_preads.len() == 0:
        return -1.0

    # Create a lookup table for raw read lengths, the key is the sequence ID.
    var lookup_raw_id = initTable[int64, int32]()
    for vals in seq_tuples_raw:
        lookup_raw_id[vals[0]] = vals[2]

    # Sum the sequence length diffs before and after error correction.
    var total_diff: int64 = 0
    for vals in seq_tuples_preads:
        # The ID pointed to by the pread name should always exist in the raw RaptorDB.
        # If it isn't there, let it fail, something went wrong.
        # Note again: this is correct: pread's NAME is the same as the rawread's ID.
        # The rawread IDs are used to name preads.
        let (pread_id, pread_name, pread_len) = vals
        let (raw_id, raw_name, raw_len) = seq_tuples_raw[parseInt(pread_name).int64]
        total_diff += (raw_len - pread_len)

    # Average value of the diffs.
    result = total_diff.float32 / seq_tuples_preads.len().float32
    return result

proc calc_length_stats*(unsorted_lengths: seq[int32]): Stats =
    if unsorted_lengths.len() == 0:
        return
    var lengths = algorithm.sorted(unsorted_lengths,
            order = algorithm.Descending)

    result.max = lengths[lengths.low.int32]
    result.min = lengths[lengths.high.int32]
    result.median = stats.percentile(lengths, 0.50)
    result.p95 = stats.percentile(lengths, 0.95)
    result.number = lengths.len().int32
    result.sum = 0
    for length in lengths:
        result.sum += length
    result.mean = 1.0 * result.sum.float32 / result.number.float32;
    result.esize = stats.esize(lengths, result.sum)
    result.nl50 = stats.len_stat(lengths, math.ceil(result.sum.float * 0.50).int32)
    result.nl95 = stats.len_stat(lengths, math.ceil(result.sum.float * 0.95).int32)

    return result

proc round2(val: float32): float =
    return parseFloat(fmt"{val:.2f}")

proc calc_stats*(genome_length, length_cutoff: int32, stats_raw, stats_seed, stats_preads: Stats, truncation: float32): string =
    var preassembled_yield: float32 = 0.0
    if stats_seed.sum > 0:
        preassembled_yield = stats_preads.sum.float / stats_seed.sum.float

    var raw_coverage: float32 = stats_raw.sum.float / genome_length.float

    var results = %*
        {
            "genome_length": genome_length,
            "length_cutoff": length_cutoff,
            "raw_reads": stats_raw.number,
            "raw_bases": stats_raw.sum,
            "raw_mean": round2(stats_raw.mean),
            "raw_n50": stats_raw.nl50[0],
            "raw_p95": stats_raw.p95,
            "raw_coverage": parseFloat(fmt"{raw_coverage:.2f}"),
            "raw_esize": round2(stats_raw.esize),
            "seed_reads": stats_seed.number,
            "seed_bases": stats_seed.sum,
            "seed_mean": round2(stats_seed.mean),
            "seed_n50": stats_seed.nl50[0],
            "seed_p95": stats_seed.p95,
            "seed_coverage": round2(stats_seed.sum.float / genome_length.float),
            "seed_esize": round2(stats_seed.esize),
            "preassembled_reads": stats_preads.number,
            "preassembled_bases": stats_preads.sum,
            "preassembled_mean": round2(stats_preads.mean),
            "preassembled_n50": stats_preads.nl50[0],
            "preassembled_p95": stats_preads.p95,
            "preassembled_coverage": round2(stats_preads.sum.float / genome_length.float),
            "preassembled_esize": round2(stats_preads.esize),
            "preassembled_yield": round2(preassembled_yield),
            # Fragmentation in HGAP4 represents the average number of preads that a seed read is broken into.
                # This is not applicable here, because we limit at most 1 fragment per seed read, so that we avoid
                # potential coverage biases and issues with chimeric reads.
            "preassembled_seed_fragmentation": 1.0,
            "preassembled_seed_truncation": round2(truncation),
        }
    return json.pretty(results, indent = 4)

proc run_from_streams*(stream_rawreads, stream_preads: Stream, genome_length, cutoff_length: int32): string =
    let tuples_raw = parse_seq_len_tuples(stream_rawreads)
    let lens_raw = get_lens_from_seq_tuples(tuples_raw)
    let stats_raw = calc_length_stats(lens_raw)

    let lens_seed = cutoff_lengths(cutoff_length, lens_raw)
    let stats_seed = calc_length_stats(lens_seed)

    let tuples_preads = parse_seq_len_tuples(stream_preads)
    let lens_preads = get_lens_from_seq_tuples(tuples_preads)
    let stats_preads = calc_length_stats(lens_preads)

    let truncation = compute_truncation(tuples_raw, tuples_preads)

    let stats = calc_stats(genome_length, cutoff_length, stats_raw, stats_seed, stats_preads, truncation)
    return stats

proc run*(rawreads_rdb_fn, preads_rdb_fn: string, genome_length, cutoff_length: int32) =
    let stream_rawreads = newFileStream(rawreads_rdb_fn, fmRead)
    defer: stream_rawreads.close()

    let stream_preads = newFileStream(preads_rdb_fn, fmRead)
    defer: stream_preads.close()

    let stats = run_from_streams(stream_rawreads, stream_preads, genome_length, cutoff_length)
    echo stats
