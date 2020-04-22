# vim: sw=4 ts=4 sts=4 tw=0 et:

## RaptorDB v0.2 specification.
## Each line is a separate entity. Order of lines does not matter. Lines start with a keyword character,
## depending on that encode:
## - Version: "V\t<float-version_num>"
## - Files: "F\t<int64-file_id>\t<str-file_path>\t<str-file_format>"
## - Sequences: "S\t<int64-seq_id>\t<str-header>\t<int64-seq_len>\t<int64-file_id>\t<int64-start_offset_in_file>\t<int64-data_length>"
## - Blocks: "B\t<int64-block_id>\t<int64-seq_id_start>\t<int64-seq_id_end>\t<int64-num_bases_in_block>"

from algorithm import nil
from math import nil
from sets import nil
from streams import nil
from strformat import fmt
from strutils import nil
from ./nuuid import nil
from ./util import nil
from random import nil
import json
import times
import sets
from streams import writeLine

type
    SequenceRecord* = object
        seq_id: int64
        header*: string
        seq_len: int64
        file_id: int64
        start_offset_in_file: int64
        data_len: int64
    FileRecord = object
        file_id*: int64
        file_path*: string
        file_format*: string
    BlockRecord = object
        block_id: int64
        seq_id_start: int64
        seq_id_end: int64
        num_bases_in_block: int64
    Db* = object
        version_major*, version_minor*: int
        files*: seq[FileRecord]
        seqs*: seq[SequenceRecord]
        blocks*: seq[BlockRecord]

    SeqLineWriter = object
        num_seq_lines: int64
        sout: File
        num_bases_written: int64
        block_size: int64
        block_start_bases: int64
        block_start_id: int64
        num_blocks: int64

proc initSeqLineWriter(sout: File, block_size_MB: int): SeqLineWriter =
    let block_size = block_size_MB * 1024 * 1024
    return SeqLineWriter(num_seq_lines: 0, sout: sout, num_bases_written: 0,
            block_size: block_size,
        block_start_bases: 0, block_start_id: 0, num_blocks: 0)
proc write_block(w: var SeqLineWriter) =
    let curr_bases_in_block = w.num_bases_written - w.block_start_bases;
    let line = fmt(
            "B\t{w.num_blocks}\t{w.block_start_id}\t{w.num_seq_lines}\t{curr_bases_in_block}\n")
    w.sout.write(line)
    inc w.num_blocks
    w.block_start_id = w.num_seq_lines
    w.block_start_bases = w.num_bases_written
proc write(w: var SeqLineWriter, split_line: seq[string]) =
    # This function updates the sequence IDs depending on the number
    # of already output sequence lines.
    var temp = split_line
    temp[1] = $w.num_seq_lines
    w.sout.write(strutils.join(temp, "\t"))
    w.sout.write('\n')
    inc w.num_seq_lines
    let nbases: int64 = strutils.parseBiggestInt(split_line[3])
    w.num_bases_written += nbases

    # Write the block if we reached the limit.
    let curr_bases_in_block = w.num_bases_written - w.block_start_bases;
    if curr_bases_in_block >= w.block_size:
        w.write_block()
proc close(w: var SeqLineWriter) =
    w.block_size = 1 # to force the last block if anything has been written
    w.write_block()

# For the blacklist filter, instead of reading the DB into memory,
# we stream.
# This might actually be slower (b/c of line-splitting into strings),
# but it will take less memory.
#
proc stream_seqs(sin, sout: File, blacklist: sets.HashSet[string]) =
    let block_size_MB = 1
    var writer = initSeqLineWriter(sout, block_size_MB)
    defer: writer.close()

    for line in lines(sin):
        if len(line) == 0:
            continue
        let keyletter = line[0]
        case keyletter
        of 'V', 'F':
            # Version and File list information.
            sout.write(line)
            sout.write('\n')
            continue
        of 'B':
            # We need to re-block because filtering will
            # reduce the size.
            continue
        of 'S':
            var fields = strutils.splitWhitespace(line)
            assert len(fields) == 7
            if not sets.contains(blacklist, fields[2]):
                writer.write(fields)
                #sout.write(line)
                #sout.write('\n')
            continue
        else:
            let msg = fmt("Unknown line '{line}'")
            #util.raiseEx(msg)
            util.log(msg)
            continue


proc filter*(blacklist_fn: string = "") =
    ## Read/write raptor-db to/from stdin/stdout.
    ## Exclude zmws in blacklist.
    # TODO: Test this. And maybe refactor with re-blocking.
    util.log("filter sans ", blacklist_fn)
    var blacklist = sets.initHashSet[string]()
    if len(blacklist_fn) > 0:
        var sin: File
        if not open(sin, blacklist_fn, fmRead):
            util.raiseEx(fmt("Failed to open blacklist file '{blacklist_fn}'"))
        defer: close(sin)
        for zmw in lines(sin):
            util.log(fmt("Skipping ({zmw})"))
            if sets.contains(blacklist, zmw):
                util.raiseEx(fmt(
                        "Found a repeat in blacklist '{blacklist_fn}': '{zmw}'\n Something is wrong!"))
            sets.incl(blacklist, zmw)
    stream_seqs(stdin, stdout, blacklist)

proc load_rdb*(sin: streams.Stream): ref Db =
    new(result)
    #var tab: char # to verify that we have read the entire "header"
    #var header: string
    var buf0: util.Headroom
    var buf1: util.Headroom

    # Delimiters should be single tabs, but we accept more in some cases.
    let v_frmt = "V %d.%d"
    let b_frmt = "B %lld %lld %lld %lld"
    let s_frmt = strutils.format("S %lld %$#s %lld %lld %lld %lld",
            (util.MAX_HEADROOM - 1))
    let f_frmt = strutils.format("F %lld %$#s %$#s",
            (util.MAX_HEADROOM - 1), (util.MAX_HEADROOM - 1))

    for line in streams.lines(sin):
        # We skip stripping, to be more strict. But we still skip totally blank lines.
        if len(line) == 0:
            continue
        elif line[0] == 'S':
            var sr: SequenceRecord
            let scanned = util.sscanf(line.cstring, s_frmt.cstring,
                addr sr.seq_id, addr buf0, addr sr.seq_len, addr sr.file_id,
                addr sr.start_offset_in_file, addr sr.data_len)
            if 6 != scanned:
                let msg = "Too few fields for '" & line & "'"
                raise newException(util.TooFewFieldsError, msg)
            util.toString(buf0, sr.header, line)
            result.seqs.add(sr)
        elif line[0] == 'B':
            var br: BlockRecord
            let scanned = util.sscanf(line.cstring, b_frmt.cstring,
                addr br.block_id, addr br.seq_id_start, addr br.seq_id_end,
                addr br.num_bases_in_block)
            if 4 != scanned:
                let msg = "Too few fields for '" & line & "'"
                raise newException(util.TooFewFieldsError, msg)
            result.blocks.add(br)
        elif line[0] == 'F':
            var fr: FileRecord
            let scanned = util.sscanf(line.cstring, f_frmt.cstring,
                addr fr.file_id, addr buf0, addr buf1)
            if 3 != scanned:
                let msg = "Too few fields for '" & line & "'"
                raise newException(util.TooFewFieldsError, msg)
            util.toString(buf0, fr.file_path, line)
            util.toString(buf1, fr.file_format, line)
            result.files.add(fr)
        elif line[0] == 'V':
            let scanned = util.sscanf(line.cstring, v_frmt.cstring,
                addr result.version_major, addr result.version_minor)
            if 2 != scanned:
                let msg = "Too few fields for '" & line & "'"
                raise newException(util.TooFewFieldsError, msg)
        else:
            let msg = "Skipping unexpected line '" & line & "'"
            util.log(msg)

proc format_version_record(version_major, version_minor: int): string =
    return "V\t{version_major}.{version_minor}".fmt
proc format_file_record(fr: FileRecord): string =
    return "F\t{fr.file_id}\t{fr.file_path}\t{fr.file_format}".fmt
proc format_seq_record(sr: SequenceRecord): string =
    return "S\t{sr.seq_id}\t{sr.header}\t{sr.seq_len}\t{sr.file_id}\t{sr.start_offset_in_file}\t{sr.data_len}".fmt
proc format_block_record(br: BlockRecord): string =
    return "B\t{br.block_id}\t{br.seq_id_start}\t{br.seq_id_end}\t{br.num_bases_in_block}".fmt

proc write_rdb*(sout: var streams.Stream; rdb: Db) =
    # If the DB is completely empty, do not write the version line.
    if (rdb.version_major, rdb.version_minor, rdb.files.len, rdb.seqs.len, rdb.blocks.len) == (0, 0, 0, 0, 0):
        return

    var line: string
    line = format_version_record(rdb.version_major, rdb.version_minor)
    sout.writeLine(line)
    for fr in rdb.files:
        line = format_file_record(fr)
        sout.writeLine(line)
    for sr in rdb.seqs:
        line = format_seq_record(sr)
        sout.writeLine(line)
    for br in rdb.blocks:
        line = format_block_record(br)
        sout.writeLine(line)

proc get_length_cutoff*(rdb_stream: streams.Stream, genome_size: int64,
        coverage: float, fail_low_cov: bool = false,
        alarms_file: string = ""): int64 =
    # The defaults here are only to help with tests.
    assert coverage > 0, fmt"Coverage needs to be > 0. Provided value: coverage = {coverage}"
    assert genome_size > 0, fmt"Genome size needs to be > 0. Provided value: genome_size = {genome_size}"

    let db = load_rdb(rdb_stream)
    let seqs = algorithm.sorted(db.seqs) do (a, b: SequenceRecord) -> int:
        return cmp(b.seq_len, a.seq_len)

    let min_desired_size = math.ceil(genome_size.float * coverage).int64
    var sum_sizes: int64 = 0
    var last_size: int64 = 0

    for sq in seqs:
        last_size = sq.seq_len
        sum_sizes += last_size
        #echo fmt"{sq.seq_len} {last_size} {sum_sizes}"
        if sum_sizes >= min_desired_size:
            return last_size # Success!

    if fail_low_cov:
        let
            min_desired_size_str = util.thousands(min_desired_size)
            sum_sizes_str = util.thousands(sum_sizes)
            msg = fmt"Not enough reads available for desired genome coverage (bases needed={min_desired_size_str} > actual={sum_sizes_str})"
        raise newException(util.GenomeCoverageError, msg)
    return last_size # This is N100, aka the min.


proc alarm(e: ref Exception, fn: string) =
    ## Write a special JSON object expected by pbcommand.models.common.
    var fout = open(fn, fmWrite)
    if nil == fout:
        util.raiseEx("Could not open '" & fn & "' for write.")
    defer: fout.close()
    let uuid = nuuid.generateUUID()
    let createdAt = $times.now().utc() #datetime.datetime.now().strftime('%Y-%m-%dT%H:%M:%S'),
                                       # this is propagated to SMRT Link UI
                                       # see PacBioAlarm class in pbcommand.models.common for details -- nat
    let alarms = %* [
        {
            "exception": $e.name,
            "info": "",
            "message": $e.msg,
            "name": $e.name,
            "severity": "ERROR",
            "owner": "TAG",
            "createdAt": createdAt,
            "id": uuid,
        }
    ]
    fout.write($alarms & "\n")


proc calc_length_cutoff*(rdb_fn: string = "rawreads.db",
        genome_size: int64 = 4600000, coverage: float = 30,
        fail_low_cov: bool = false, alarms_fn: string = "") =
    ## Perform a linear pass on an overlap file, and determine rough clipping coordinates to 'correct' reads.
    ## Write integer to stdout.
    try:
        let content = system.readFile(rdb_fn) # read into memory b/c streams lib is slow
        var sin = streams.newStringStream(content)
        defer: streams.close(sin)

        let cutoff = get_length_cutoff(sin, genome_size, coverage, fail_low_cov, alarms_fn)
        stdout.write($cutoff & "\n")
    except Exception as exc:
        if "" != alarms_fn:
            alarm(exc, alarms_fn)
        raise


proc rdb_reblock(rdb: var Db, block_size: int64) =
    # Clear out the blocks.
    rdb.blocks = @[]

    # Empty input, do nothing.
    if rdb.seqs.len() == 0:
        return

    # Initialize the buffer for a new block.
    var new_block: BlockRecord
    new_block.block_id = rdb.blocks.len
    new_block.seq_id_start = 0
    new_block.seq_id_end = 0
    new_block.num_bases_in_block = 0

    # Re-enumerate the sequence IDs, and assign them to blocks.
    for i in 0..<rdb.seqs.len:
        rdb.seqs[i].seq_id = i

        # Extend the block.
        new_block.seq_id_end = i + 1
        new_block.num_bases_in_block += rdb.seqs[i].seq_len

        # Start a new block if needed.
        if block_size >= 0 and new_block.num_bases_in_block >= block_size:
            rdb.blocks.add(new_block)
            new_block = BlockRecord()

            new_block.block_id = rdb.blocks.len
            new_block.seq_id_start = i + 1
            new_block.seq_id_end = i + 1
            new_block.num_bases_in_block = 0

    # Add the last block.
    if new_block.seq_id_end > new_block.seq_id_start:
        rdb.blocks.add(new_block)

proc split_movie_name(header: string): tuple[movie_name: string, zmw_id: string, zmw_pos: string] =
    let sp_header = strutils.split(header, '/')
    var ret = (sp_header[0], sp_header[1], sp_header[2])
    return ret

proc get_subsampled_rdb*(rdb: ref Db, genome_size: int64,
        coverage: float, use_umc: bool, random_seed: int64, block_size: int64, fail_low_cov: bool = false,
        alarms_file: string = ""): Db =
    assert coverage > 0, fmt"Coverage needs to be > 0. Provided value: coverage = {coverage}"
    assert genome_size > 0, fmt"Genome size needs to be > 0. Provided value: genome_size = {genome_size}"

    # Reproducibility, if a valid seed is given.
    if random_seed > 0:
        random.randomize(random_seed)

    let min_desired_size: int64 = int64(float(genome_size) * coverage)

    # Prepare an output DB and copy info which will be unchanged.
    var ret_rdb: Db
    ret_rdb.files = rdb.files
    ret_rdb.version_major = rdb.version_major
    ret_rdb.version_minor = rdb.version_minor

    var selected_set = initHashSet[string]()
    var sum_sizes: int64 = 0

    # Create a random permutation.
    var permuted_seqs = rdb.seqs
    random.shuffle(permuted_seqs)

    # Select the names
    for seq in permuted_seqs:
        var name = seq.header

        # If use_umc is true, we will output all subreads for a particular selected ZMW.
        if use_umc == true:
            let (movie_name, zmw_id, zmw_pos) = split_movie_name(seq.header)
            name = "{movie_name}/{zmw_id}".fmt

        # Skip duplicates. Important for ZMW subsampling.
        if selected_set.contains(name):
            continue
        # Accumulate the new sequence.
        selected_set.incl(name)
        sum_sizes += seq.seq_len
        if sum_sizes >= min_desired_size:
            break

    # Check if we should raise.
    if fail_low_cov and sum_sizes < min_desired_size:
        let
            min_desired_size_str = util.thousands(min_desired_size)
            sum_sizes_str = util.thousands(sum_sizes)
            msg = fmt"Not enough reads available to downsample to desired genome coverage (bases needed={min_desired_size_str} > actual={sum_sizes_str})"
        raise newException(util.GenomeCoverageError, msg)

    # Second pass, pick the actual sequences.
    for seq in rdb.seqs:
        var name = seq.header
        if use_umc == true:
            let (movie_name, zmw_id, zmw_pos) = split_movie_name(seq.header)
            name = "{movie_name}/{zmw_id}".fmt
        if name in selected_set:
            ret_rdb.seqs.add(seq)

    # Re-enumerate the sequences and form the blocks, in-place.
    rdb_reblock(ret_rdb, block_size)

    return ret_rdb

proc run_subsample*(rdb_fn: string, genome_size: int64, coverage: float, use_umc: bool, random_seed: int64, block_size_mb: float, alarms_fn: string = "") =
    let block_size_in_bytes = int64(block_size_mb * 1024.0 * 1024.0)

    try:
        let content = system.readFile(rdb_fn) # read into memory b/c streams lib is slow
        var sin = streams.newStringStream(content)
        defer: streams.close(sin)

        var sout: streams.Stream = streams.newFileStream(stdout)
        defer: streams.close(sout)

        let rdb = load_rdb(sin)
        let subsampled_rdb = get_subsampled_rdb(rdb, genome_size, coverage, use_umc, random_seed, block_size_in_bytes, false, alarms_fn)
        write_rdb(sout, subsampled_rdb)

    except Exception as exc:
        if "" != alarms_fn:
            alarm(exc, alarms_fn)
        raise

