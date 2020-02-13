# vim: sw=4 ts=4 sts=4 tw=0 et:
from algorithm import nil
import streams
import strformat
from json import `%*`
import overlapParser
import std/stats
import system
import os
from strutils import toLowerAscii

proc median*(xs: seq[int]): float =
    if xs.len() == 0: return 0.0
    var ys = xs
    algorithm.sort(ys, system.cmp[int])
    # The seq.high is the index of the last element in the seq, i.e. one less than len.
    return 0.5 * float(ys[ys.high div 2] + ys[ys.len div 2]) 

proc count_overlaps*(pile: seq[overlapParser.Overlap]): (int, int) =
    var count5, count3: int = 0
    for ov in pile:
        if ov.Astart == 0:
            inc(count5)
        if ov.Aend == ov.Alen:
            inc(count3)
    return (count5, count3)

proc collect_from_stream*(inStream: Stream): seq[int] =
    # Collect both 5' and 3' stats because coverage may vary.
    var allCounts: seq[int]
    for pile in overlapParser.getNextPile(inStream):
        let (count5, count3) = count_overlaps(pile)
        allCounts.add(count5)
        allCounts.add(count3)
    return allCounts

proc collect_from_fofn*(in_fn_list: seq[string]): seq[int] =
    var allCounts: seq[int]
    for fn in in_fn_list:
        stderr.writeLine("Processing: {fn}".fmt)
        let inStream = newFileStream(fn, fmRead)
        defer: inStream.close()
        let counts = collect_from_stream(inStream)
        allCounts.add(counts)
    return allCounts

proc calc_stats*(counts: seq[int]): tuple[mean, stddev, median: float, min, max, n: int] =
    var rs: RunningStat
    rs.push(counts)
    let median_val = median(counts)
    return (mean: rs.mean(), stddev: rs.standardDeviation(), median: median_val, min: int(rs.min), max: int(rs.max), n: rs.n)

proc parse_fofn*(in_fn: string): seq[string] =
    let inStream = newFileStream(in_fn, fmRead)
    defer: inStream.close()
    var line: string
    var ret: seq[string]
    while(readLine(inStream, line)):
        ret.add(line)
    return ret

proc run*(in_fn: string) =
    let fileSplit = os.splitFile(in_fn)
    let is_fofn = (toLowerAscii(fileSplit.ext) == ".fofn")

    var in_fn_list = @[in_fn]
    if is_fofn: in_fn_list = parse_fofn(in_fn)

    let counts = collect_from_fofn(in_fn_list)
    let (mean, sd, med, min, max, n) = calc_stats(counts)

    var results = %*
        {
            "mean": mean,
            "std": sd,
            "median": med,
            "min": min,
            "max": max,
            "n": n,
        }
    echo json.pretty(results, indent = 4)
