# vim: sw=4 ts=4 sts=4 tw=0 et:
from ./util import log
from algorithm import nil
from os import nil
from math import nil
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

proc countLines*(sin: streams.Stream): int =
    # Simply count the number of lines in this file.
    # (Each line is a read.)
    var line: string
    while streams.readLine(sin, line):
        result += 1

proc getBlockIdFromFilename*(fn: string, ext = ".reads"): int =
    # foo.99.reads => 99
    # This could use re, but to avoid requiring pcre we
    # simply search for ".reads".
    let lastDot = len(fn) - len(ext)
    let next2LastDot = strutils.rfind(fn, '.', 0, lastDot-1)
    assert next2LastDot != -1, fn
    #echo "next:", next2LastDot, " ", lastDot, " '", fn[next2LastDot+1 ..< lastDot], "'"
    return strutils.parseInt(fn[next2LastDot+1 ..< lastDot])

proc countLines(prefix, ext: string): BlockWeights =
    # Find all files called {prefix}(.*){suffix} in the filesystem,
    # and count the number of reads in each.
    # E.g. "blocks.*.m4"
    # In the result, the index is the block-id.
    #let counts = tables.toCountTable(sequtils.toSeq(tables.values(r2c)))
    #var
    #    contigs = sequtils.toSeq(tables.keys(counts))
    #    sizes: seq[int]
    #algorithm.sort(contigs)
    let pattern = prefix & "*" & ext
    var counts = tables.initCountTable[int]()
    for fn in os.walkFiles(pattern):
        let sin = streams.openFileStream(fn)
        let count = countLines(sin)
        streams.close(sin)
        log("Found {count} lines in '{fn}'.".fmt)
        let block_id = getBlockIdFromFilename(fn, ext)
        tables.inc(counts, block_id, count)
    for block_id in 0 ..< len(counts): # TODO: What if input chunks are not consecutive?
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
        fout.close()
    return shards

proc partitionBlocks(prefix: string, blockWeights: BlockWeights, nShards: int): seq[seq[int]] =
    # Combine blocks into at most nShards subsets, weighted by the
    # number of reads for each block (which map to one or more contigs).
    # Assume block_ids are 0 .. len(blockWeights)-1.
    # This version allows non-consecutive block_ids for each shard.
    # len(result) will be <= nShards

    let shards: seq[seq[int]] = util.partitionWeighted(nShards, blockWeights)
    var shard_sizes: seq[int]
    for s in shards:
        shard_sizes.add(s.len())
    log("shard sizes:{shard_sizes}".fmt)
    for shard_id in 0 ..< shards.len():
        let fn = prefix & "." & $shard_id & ".block_ids"
        log("shard:{fn}".fmt)
        var fout = open(fn, fmWrite)
        for block_id in shards[shard_id]:
            fout.writeLine(block_id)
        fout.close()
    return shards

type
    PancakeRange* = object
        t*, qs*, qe*: int # target, querystart, queryend

proc size*(a: PancakeRange): int =
    return a.qe - a.qs

proc splitPancakeRange*(a: PancakeRange, leftSize: int): tuple[left: PancakeRange, right: PancakeRange] =
    # Split into 2 shorter ranges covering the whole.
    # Left is at most leftSize; right is the rest.
    assert leftSize <= a.size()
    let
        left = PancakeRange(t: a.t, qs: a.qs, qe: (a.qs + leftSize))
        right = PancakeRange(t: a.t, qs: left.qe, qe: a.qe)
    assert left.size() == leftSize
    assert left.size() + right.size() == a.size()
    return (left, right)

proc splitPancakeRangeRight*(a: PancakeRange, rightSize: int): tuple[left: PancakeRange, right: PancakeRange] =
    # Split into 2 shorter ranges covering the whole.
    # Right is at most rightSize; left is the rest.
    return splitPancakeRange(a, a.size - rightSize)

proc `$`*(a: PancakeRange): string =
    return "{a.t} {a.qs} {a.qe}".fmt

proc shardMatrix*(nrows, ncols, nShards: int): seq[seq[PancakeRange]] =
    ## Note order of arguments.
    # Initialize stack w/ n complete rows.
    var stack = newSeq[PancakeRange](nrows)
    var summed = 0
    for i in 0 ..< nrows:
        let pr = PancakeRange(t: i, qs: 0, qe: ncols)
        summed += pr.size()
        stack[nrows - i - 1] = pr
    # For each shard, consume total/nshards comparisons, splitting
    # as necessary.
    var
        total = nrows*ncols
        remaining = total
        needed = 0
        i = 0
    assert total == summed
    while remaining > 0:
        needed = math.ceil(remaining/(nShards - result.len())).int # rounded up
        while needed > 0:
            result.setLen(i + 1)
            #log("i:{i} rem:{remaining} needed:{needed} max:{nShards}".fmt)
            var pr = stack.pop()
            if pr.size() > needed:
                # Snake back and forth for better data-locality.
                if (pr.t and 1) == 0:
                    # even row
                    let (left, right) = pr.splitPancakeRange(needed)
                    stack.add(right) # Push back what we did not yet need.
                    pr = left
                else:
                    let (left, right) = pr.splitPancakeRangeRight(needed)
                    stack.add(left) # Push back what we did not yet need.
                    pr = right
                    # odd row
            result[i].add(pr)
            #log(" result[i]:{result[i]}".fmt)
            remaining -= pr.size()
            needed -= pr.size()
        i += 1
    #echo "result=" & $result
    assert len(result) <= nShards
    return result

proc shardMatrixColumns*(nrows, ncols, nShards: int): seq[seq[PancakeRange]] =
    ## Note order of arguments.
    if nrows == 0 or ncols == 0:
        return
    # Initialize stack w/ a single, complete row.
    # (This lets us re-use some logic.)
    # (We will duplicate rows later.)
    let
        pr = PancakeRange(t: 0, qs: 0, qe: ncols)
        summed = pr.size()
    var stack = @[pr]

    # For each shard, consume total/nshards comparisons, splitting
    # as necessary.
    var
        total = 1*ncols
        remaining = total
        needed = 0
        i = 0
    assert total == summed
    while remaining > 0:
        needed = math.ceil(remaining/(nShards - result.len())).int # rounded up
        while needed > 0:
            result.setLen(i + 1)
            #log("i:{i} rem:{remaining} needed:{needed} max:{nShards}".fmt)
            var pr = stack.pop()
            if pr.size() > needed:
                let (left, right) = pr.splitPancakeRange(needed)
                stack.add(right) # Push back what we did not yet need.
                pr = left
            let
                qs = pr.qs
                qe = pr.qe
            for r in 0 ..< nrows:
                let row = PancakeRange(t: r, qs: qs, qe: qe)
                result[i].add(row)
            #log(" result[i]:{result[i]}".fmt)
            remaining -= pr.size()
            needed -= pr.size()
        i += 1
    #echo "result=" & $result
    assert len(result) <= nShards
    return result

proc shardUpperTriangular*(n: int, nShards: int): seq[seq[PancakeRange]] =
    # Combine upper-triangular ranges of blocks into at most nShards subsets.
    # We need 0 vs. [0..n), 1 vs. [1..n), 2 vs. [2..n), etc.
    # Think of all those comparisons (there are C==n*(n+1)/2) in a straight line,
    # and take roughly C/nShards for each shard. Then turn those back into
    # contiguous rows for efficiency in Pancake.
    # (See tests for examples.)
    # len(result) will be <= nShards

    var stack = newSeq[PancakeRange](n)
    # Initialize stack w/ n complete rows.
    var summed = 0
    for i in 0 ..< n:
        let pr = PancakeRange(t: i, qs: i, qe: n)
        summed += pr.size()
        stack[n - i - 1] = pr
    # For each shard, consume total/nshards comparisons, splitting
    # as necessary.
    var
        total = (n * (n+1) / 2).int
        remaining = total
        needed = 0
        i = 0
    assert total == summed
    while remaining > 0:
        needed = math.ceil(remaining/(nShards - result.len())).int # rounded up
        while needed > 0:
            result.setLen(i + 1)
            #log("i:{i} rem:{remaining} needed:{needed} max:{nShards}".fmt)
            var pr = stack.pop()
            if pr.size() > needed:
                let (left, right) = pr.splitPancakeRange(needed)
                stack.add(right) # Push back what we did not yet need.
                pr = left
            result[i].add(pr)
            #log(" result[i]:{result[i]}".fmt)
            remaining -= pr.size()
            needed -= pr.size()
        i += 1
    #echo "result=" & $result
    assert len(result) <= nShards
    return result

proc combineShards(prefix: string, shards: seq[seq[PancakeRange]]) =
    var shard_sizes: seq[int]
    for s in shards:
        shard_sizes.add(s.len())
    log("shard sizes:{shard_sizes}".fmt)
    for shard_id in 0 ..< shards.len():
        let fn = prefix & "." & $shard_id & ".pancake_ranges"
        log("shard:{fn}".fmt)
        var fout = open(fn, fmWrite)
        for block_range in shards[shard_id]:
            fout.writeLine(block_range)
        fout.close()

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

proc split_paf*(max_nshards: int, shard_prefix = "shard", block_prefix = "block",
        in_read_to_contig_fn = "sorted.read_to_contig.csv", out_ids_fn = "all_shard_ids",
        mb_per_block: int,
        blacklist_fn: string = "",
        in_fai_fns: seq[string]) =
    ## The trailing list of fasta.fai filenames are FASTA index files.
    ## They will be used to split the shards somewhat evenly.
    ## (Used to shard the polishing jobs.)

    var blacklist = initHashSet[string]()
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
    let shards = combineBlocks(shard_prefix, countLines(block_prefix&".", ".reads"), max_nshards)
    if out_ids_fn != "":
        var fout = open(out_ids_fn, fmWrite)
        for shard_id in 0 ..< len(shards):
            fout.writeLine(shard_id)
        fout.close()

proc split*(max_nshards: int, shard_prefix = "shard", block_prefix = "block",
        in_read_to_contig_fn = "sorted.read_to_contig.csv", out_ids_fn = "all_shard_ids",
        mb_per_block: int,
        blacklist_fn: string = "",
        in_fai_fns: seq[string]) =
    ## The trailing list of fasta.fai filenames are FASTA index files.
    ## They will be used to split the shards somewhat evenly.
    ## (Used to shard the polishing jobs.)

    var blacklist = initHashSet[string]()
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
    let shards = combineBlocks(shard_prefix, countLines(block_prefix&".", ".reads"), max_nshards)
    if out_ids_fn != "":
        var fout = open(out_ids_fn, fmWrite)
        for shard_id in 0 ..< len(shards):
            fout.writeLine(shard_id)
        fout.close()

proc shard_blocks_m4*(max_nshards: int, shard_prefix = "shard", block_prefix = "block",
        out_ids_fn = "all_shard_ids") =
    ## Given several {block_prefix}.(block_id).m4 files,
    ## create up to {max_nshards} files that each contain a list
    ## of block_ids (one per line).
    ## For now, they are balanced by the number of reads in each .m4 file.
    ## (Later, the contents of each shard will be processed linearly, one block at a time,
    ## on a given compute node.)
    ## (Used to shard the phasing jobs.)
    let shards = partitionBlocks(shard_prefix, countLines(block_prefix&".", ".m4"), max_nshards)
    if out_ids_fn != "":
        var fout = open(out_ids_fn, fmWrite)
        for shard_id in 0 ..< len(shards):
            fout.writeLine(shard_id)
        fout.close()

proc shard_ovl_asym*(max_nshards: int, shard_prefix = "shard", n: int, out_ids_fn = "") =
    ## (Used to shard the asymmetric overlap jobs.)
    log("shard_ovl_asym: max_nshards={max_nshards} n={n}".fmt)
    let shards: seq[seq[PancakeRange]] = shardUpperTriangular(n = n, nShards = max_nshards)
    combineShards(prefix = shard_prefix, shards)
    if out_ids_fn != "":
        var fout = open(out_ids_fn, fmWrite)
        for shard_id in 0 ..< len(shards):
            fout.writeLine(shard_id)
        fout.close()

proc shard_mapping*(max_nshards: int, shard_prefix = "shard", n_query_blocks, n_target_blocks: int, out_ids_fn = "") =
    ## Generate comparisons for nq-by-nt matrix.
    ## (Used to shard the purge_dups overlap jobs, contigs vs. reads.)
    log("shard_mapping: max_nshards={max_nshards} n_query={n_query_blocks} n_target={n_target_blocks}".fmt)
    let shards = shardMatrixColumns(nrows = n_target_blocks, ncols = n_query_blocks, nShards = max_nshards)
    combineShards(prefix = shard_prefix, shards)
    if out_ids_fn != "":
        var fout = open(out_ids_fn, fmWrite)
        for shard_id in 0 ..< len(shards):
            fout.writeLine(shard_id)
        fout.close()

proc prepare*(max_nshards: int, shard_prefix = "shard_id", block_prefix = "block_id", out_ids_fn = "") =
    ## DEPRECATED
    log("prepare {max_nshards} shard:{shard_prefix} block:{block_prefix}".fmt)
    let shards = combineBlocks(shard_prefix, countLines(block_prefix&".", ".reads"), max_nshards)
    if out_ids_fn != "":
        var fout = open(out_ids_fn, fmWrite)
        for shard_id in 0 ..< len(shards):
            fout.writeLine(shard_id)
        fout.close()
