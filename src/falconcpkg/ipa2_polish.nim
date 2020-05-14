# vim: sw=4 ts=4 sts=4 tw=0 et:
from ./util import log
from algorithm import nil
from os import nil
from parsecsv import nil
from sequtils import nil
from streams import nil
from strformat import fmt
from strutils import nil
import hts
import tables
import sets

type
    BlockWeights = seq[int]

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

proc countBlocks(prefix: string): BlockWeights =
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

proc combineBlocks(prefix: string, blockWeights: BlockWeights, nShards: int): seq[int] =
    # Combine blocks into at most nShards subsets, weighted by the
    # number of reads for each block (which map to one or more contigs).
    # Assume block_ids are 0 .. len(blockWeights)-1.
    # len(result) will be <= nShards

    let shards: seq[int] = util.splitWeighted(nShards, blockWeights)
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

type
    Contig2Reads = tables.TableRef[string, seq[string]]

proc get_Contig2Reads(sin: streams.Stream, in_r2c_fn: string, contig2len: tables.TableRef[string, int]): Contig2Reads =
    result = tables.newTable[string, seq[string]]()
    var parser: parsecsv.CsvParser
    parsecsv.open(parser, sin, filename = in_r2c_fn, separator = ' ', skipInitialSpace = true)
    while parsecsv.readRow(parser, 2):
        if not contig2len.haskey(parser.row[1]):
            continue
        tables.mgetOrPut(result, parser.row[1], @[]).add(parser.row[0])

proc combineContigs(target_mb: int, contigs: seq[string], contig2len: tables.TableRef[string, int]): seq[seq[string]] =
    # Combine contigs into subsets, where each has at least mb MegaBases
    # of contigs.
    if len(contigs) != tables.len(contig2len):
        let msg = "len(contigs)={len(contigs)} != len(contig2len)={len(contig2len)}".fmt
        util.raiseEx(msg)
    var weights: seq[int64]
    for contig_idx in 0 ..< len(contigs):
        weights.add(contig2len[contigs[contig_idx]])
    let contig_indices: seq[seq[int]] = util.combineToTarget(target_mb * 1_000_000, weights)
    for block_idx in 0 ..< len(contig_indices):
        assert 0 != len(contig_indices[block_idx])
        result.add(@[])
        for contig_idx in contig_indices[block_idx]:
            result[block_idx].add(contigs[contig_idx])

proc writeBlocksMulti(contig2len: tables.TableRef[string, int], contig2reads: Contig2Reads, prefix: string, mb: int) =
    # We want multiple contigs per block, as many as slightly over mb
    # MegaBases of contig lengths.
    # For now, all contigs are in *single* block.
    var contigs = sequtils.toSeq(tables.keys(contig2reads))
    algorithm.sort(contigs) # Sort only for stable tests. Random order is fine too.

    let blocks = combineContigs(mb, contigs, contig2len)

    var count = 0
    for bloke in blocks:
        let
            ctg_id_fn = "{prefix}.{count}.ctg_id".fmt
            reads_fn = "{prefix}.{count}.reads".fmt
        var
            ctg_id_fout = system.open(ctg_id_fn, fmWrite)
            reads_fout = system.open(reads_fn, fmWrite)
        for contig in bloke:
            ctg_id_fout.writeLine(contig)
            let
                reads = contig2reads[contig]
            for read in reads:
                reads_fout.writeLine(read)
        system.close(reads_fout)
        system.close(ctg_id_fout)
        count += 1

proc writeBlocksLikeAwk(contig2reads: Contig2Reads, prefix: string) =
    discard """
        awk -v output_blocks="all_block_ids" 'BEGIN { prev = ""; count = 0; fn_blocks = output_blocks} { if ($2 != prev) { fn_reads = "block."count".reads"; fn_ctg =  "block."count".ctg_id"; print $2 > fn_ctg; print count > fn_blocks; count++} {print $1 > fn_reads} prev = $2; }' sorted.read_to_contig.csv
    """
    var contigs = sequtils.toSeq(tables.keys(contig2reads))
    algorithm.sort(contigs) # Sort only for stable tests.
    for count in 0 ..< len(contigs):
        let
            contig = contigs[count]
            reads = contig2reads[contig]
            #count = reads.len()
            ctg_id_fn = "block.{count}.ctg_id".fmt
            reads_fn = "block.{count}.reads".fmt
        block:
            var fout = system.open(ctg_id_fn, fmWrite)
            fout.writeLine(contig)
            system.close(fout)
        var fout = system.open(reads_fn, fmWrite)
        for read in reads:
            fout.writeLine(read)
        system.close(fout)

proc updateChromLens(chrom2len: tables.TableRef[string, int], fai_fn: string, blacklist: HashSet[string]) =
    # Raise if a new chrom already exists in this table.
    var fin: hts.Fai
    let fn = fai_fn[0 ..< ^4]
    if not open(fin, fn):
        util.raiseEx("Could not open '{fai_fn}' ({fn})".fmt)
    let nchroms = hts.len(fin)
    log("{nchroms} reads in '{fn}'".fmt)
    for i in 0 ..< nchroms:
        let chrom: string = fin[i]
        # echo " chrom:", chrom
        if chrom in chrom2len:
            let
                n: int = chrom2len[chrom]
                msg = "Error: chrom '{chrom}' was already seen {n}. Read names in fasta input files must be unique!".fmt
            util.raiseEx(msg)
        let n = hts.chrom_len(fin, chrom)
        if blacklist.contains(chrom):
            log("Skipping blacklisted sequence: chrom = {chrom}, len = {n}".fmt)
            continue
        chrom2len[chrom] = n
    #hts.destroy_fai(fin) # We may need to activate destructors to do this properly. Oh, well.

proc loadSet(sin: streams.Stream): HashSet[string] =
    var line: string
    while streams.readLine(sin, line):
        result.incl(line)

proc split*(max_nshards: int, shard_prefix = "shard", block_prefix = "block",
        in_read_to_contig_fn = "sorted.read_to_contig.csv", out_ids_fn = "all_shard_ids",
        mb_per_block: int,
        blacklist_fn: string = "",
        in_fai_fns: seq[string]) =
    ## The trailing list of fasta.fai filenames are FASTA index files.
    ## They will be used to split the shards somewhat evenly.

    var blacklist = initSet[string]()
    if len(blacklist_fn) != 0:
        log("Loading the blacklist from '{blacklist_fn}'.".fmt)
        var sin = streams.openFileStream(blacklist_fn)
        blacklist = loadSet(sin)
        streams.close(sin)
    log("Blacklist contains {len(blacklist)} elements.".fmt)

    log("split {max_nshards} shard:{shard_prefix} block:{block_prefix} in:'{in_read_to_contig_fn}' out:'{out_ids_fn}'".fmt)
    var chrom2len = tables.newTable[string, int]()
    for fn in in_fai_fns:
        chrom2len.updateChromLens(fn, blacklist)
    var sin = streams.openFileStream(in_read_to_contig_fn)
    let contig2reads = get_Contig2Reads(sin, in_read_to_contig_fn, chrom2len)
    streams.close(sin)
    writeBlocksMulti(chrom2len, contig2reads, block_prefix, mb_per_block)
    let shards = combineBlocks(shard_prefix, countBlocks(block_prefix), max_nshards)
    if out_ids_fn != "":
        var fout = open(out_ids_fn, fmWrite)
        for shard_id in 0 ..< len(shards):
            fout.writeLine(shard_id)
        fout.close()

proc prepare*(max_nshards: int, shard_prefix = "shard_id", block_prefix = "block_id", out_ids_fn = "") =
    log("prepare {max_nshards} shard:{shard_prefix} block:{block_prefix}".fmt)
    let shards = combineBlocks(shard_prefix, countBlocks(block_prefix), max_nshards)
    if out_ids_fn != "":
        var fout = open(out_ids_fn, fmWrite)
        for shard_id in 0 ..< len(shards):
            fout.writeLine(shard_id)
        fout.close()
