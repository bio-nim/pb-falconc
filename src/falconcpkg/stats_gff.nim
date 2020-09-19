# vim: sw=4 ts=4 sts=4 tw=0 et:
import json
import streams
import strformat
import gff_parser
import tables

proc summarize_gff_coverage*(sin: Stream): Table[string, Table[string, string]] =
    type
        SeqCounts = object
            sum_cov, sum_span: float64
            cov_min, cov_max: int
            cov_mean: float64
    var seq_lookup = initTable[string, SeqCounts]()

    # Compute the mean over all windows.
    for record in yield_gff_record(sin):
        if not (record.seq_name in seq_lookup):
            seq_lookup[record.seq_name] = SeqCounts(sum_cov: 0.0, sum_span: 0.0, cov_min: high(int), cov_max: -1, cov_mean: 0)
        var bucket = addr seq_lookup[record.seq_name]
        bucket.sum_cov += float64(record.seq_end - record.seq_start) * float64(record.cov_mean)
        bucket.sum_span += float64(record.seq_end - record.seq_start)
        bucket.cov_mean = bucket.sum_cov / bucket.sum_span
        bucket.cov_min = min(bucket.cov_min, record.cov_min)
        bucket.cov_max = max(bucket.cov_max, record.cov_max)

    # Format the results dict.
    var results = initTable[string, Table[string, string]]()
    for key, val in seq_lookup:
        results[key] = initTable[string, string]()
        results[key]["cov_min"] = "{val.cov_min}".fmt
        results[key]["cov_max"] = "{val.cov_max}".fmt
        results[key]["cov_mean"] = "{val.cov_mean}".fmt

    return results

proc run*(gff_fn: string) =
    let gff_stream = newFileStream(gff_fn, fmRead)
    defer: gff_stream.close()

    let results = summarize_gff_coverage(gff_stream)
    echo json.pretty(parseJson($results), indent = 4)

