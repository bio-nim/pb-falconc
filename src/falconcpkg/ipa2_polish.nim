# vim: sw=4 ts=4 sts=4 tw=0 et:
from ./util import log
#from algorithm import nil
from os import nil
#from sequtils import nil
from streams import nil
from strformat import fmt
from strutils import nil
import tables

type
    BlockCounts = seq[int]

proc countReadsInBlock*(sin: streams.Stream): int =
    # Simply count the number of lines in this file.
    # (Each line is a read.)
    var line: string
    while streams.readLine(sin, line):
        result += 1

proc getBlockIdFromFilename*(fn: string): int =
    # foo.99.reads => 99
    # This could use re, but to avoid requiring pcre we
    # simply search for ".reads".
    let lastDot = len(fn) - 6
    let next2LastDot = strutils.rfind(fn, '.', 0, lastDot-1)
    assert next2LastDot != -1, fn
    #echo "next:", next2LastDot, " ", lastDot, " '", fn[next2LastDot+1 ..< lastDot], "'"
    return strutils.parseInt(fn[next2LastDot+1 ..< lastDot])

proc countBlocks(prefix: string): BlockCounts =
    # Find all files called prefix.*.reads in the filesystem,
    # and count the number of reads in each.
    # In the result, the index is the block-id.
    #let counts = tables.toCountTable(sequtils.toSeq(tables.values(r2c)))
    #var
    #    contigs = sequtils.toSeq(tables.keys(counts))
    #    sizes: seq[int]
    #algorithm.sort(contigs)
    let pattern = prefix & ".*.reads"
    var counts = tables.initCountTable[int]()
    for fn in os.walkFiles(pattern):
        let count = countReadsInBlock(streams.openFileStream(fn))
        log("Found {count} reads in {fn}.".fmt)
        let block_id = getBlockIdFromFilename(fn)
        tables.inc(counts, block_id, count)
    for block_id in 0 ..< len(counts):
        result.add(counts[block_id])

proc shardBlocks(prefix: string, block_counts: BlockCounts, nShards: int): seq[int] =
    # Shard contigs into at most nShards subsets, weighted by the
    # number of reads for each block (which is a contig).
    # Assume block_ids are 0 .. len(block_counts)-1.
    # len(result) will be <= nShards

    let shards: seq[int] = util.splitWeighted(nShards, blockCounts)
    log("shards:{shards}".fmt)
    var block_id = 0
    for shard_id in 0 ..< len(shards):
        let fn = prefix & "." & $shard_id & ".block_ids"
        log("shard:{fn}".fmt)
        var fout = open(fn, fmWrite)
        var count = shards[shard_id]
        while count > 0:
            fout.writeLine(block_id)
            block_id += 1
            count -= 1
    return shards

proc split*(max_nshards: int, shard_prefix = "shard_id", block_prefix = "block_id", out_ids_fn = "") =
    log("prepare {max_nshards} shard:{shard_prefix} block:{block_prefix}".fmt)
    let shards = shardBlocks(shard_prefix, countBlocks(block_prefix), max_nshards)
    if out_ids_fn != "":
        var fout = open(out_ids_fn, fmWrite)
        for shard_id in 0 ..< len(shards):
            fout.writeLine(shard_id)
        fout.close()
