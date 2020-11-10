# vim: sw=4 ts=4 sts=4 tw=0 et:

import falconcpkg/ipa2_construct_config
import falconcpkg/ipa2_polish
import falconcpkg/ipa2Header2Falcon
import unittest, streams, tables

let expected_default_bash = """
config_autocomp_max_cov=1
config_block_size=4096
config_cleanup=0
config_coverage=0
config_erc_min_idt=99.9
config_existing_db_prefix=''
config_genome_size=0
config_m4filt_high_copy_sample_rate=1.0
config_max_polish_block_mb=100
config_ovl_filter_opt='--max-diff 80 --max-cov 100 --min-cov 1 --bestn 10 --min-len 4000 --gapFilt --minDepth 4'
config_ovl_min_idt=98
config_ovl_min_len=1000
config_ovl_opt=''
config_phase_run=1
config_phasing_opt=''
config_phasing_split_opt='--split-type noverlaps --limit 3000000'
config_polish_min_len=50000
config_polish_run=1
config_purge_dups_calcuts=''
config_purge_dups_get_seqs=''
config_purge_dups_run=0
config_purge_map_opt='--min-map-len 1000 --min-idt 98.0 --bestn 5'
config_seeddb_opt='-k 32 -w 80 --space 2'
config_seqdb_opt='--compression 0'
config_use_hpc=0
config_use_seq_ids=0
"""

let expected_default_json = """
{
    "config_autocomp_max_cov": "1",
    "config_block_size": "4096",
    "config_cleanup": "0",
    "config_coverage": "0",
    "config_erc_min_idt": "99.9",
    "config_existing_db_prefix": "",
    "config_genome_size": "0",
    "config_m4filt_high_copy_sample_rate": "1.0",
    "config_max_polish_block_mb": "100",
    "config_ovl_filter_opt": "--max-diff 80 --max-cov 100 --min-cov 1 --bestn 10 --min-len 4000 --gapFilt --minDepth 4",
    "config_ovl_min_idt": "98",
    "config_ovl_min_len": "1000",
    "config_ovl_opt": "",
    "config_phase_run": "1",
    "config_phasing_opt": "",
    "config_phasing_split_opt": "--split-type noverlaps --limit 3000000",
    "config_polish_min_len": "50000",
    "config_polish_run": "1",
    "config_purge_dups_calcuts": "",
    "config_purge_dups_get_seqs": "",
    "config_purge_dups_run": "0",
    "config_purge_map_opt": "--min-map-len 1000 --min-idt 98.0 --bestn 5",
    "config_seeddb_opt": "-k 32 -w 80 --space 2",
    "config_seqdb_opt": "--compression 0",
    "config_use_hpc": "0",
    "config_use_seq_ids": "0"
}
"""

suite "ipa2_construct_config":
    test "write bash":
        var
            outs = streams.newStringStream()
            ins = streams.newStringStream()
        run(outs, ins, "", 'b', true)
        outs.setPosition(0)
        check outs.readAll() == expected_default_bash

    test "write json":
        var
            outs = streams.newStringStream()
            ins = streams.newStringStream()
        run(outs, ins, "", 'j', true)
        outs.setPosition(0)
        check outs.readAll() == expected_default_json

    test "unset":
        let cfg = parse("config_coverage=")
        check cfg["config_coverage"] == ""
        check cfg["config_block_size"] == "4096"

    test "override":
        let cfg = parse("config_coverage=11")
        check cfg["config_coverage"] == "11"
        check cfg["config_block_size"] == "4096"

    test "separated by semicolon":
        let cfg = parse("config_coverage=11; config_genome_size=12")
        check cfg["config_coverage"] == "11"
        check cfg["config_genome_size"] == "12"
        check cfg["config_block_size"] == "4096"

    test "separated by newline":
        let cfg = parse("config_coverage=11\nconfig_genome_size=12")
        check cfg["config_coverage"] == "11"
        check cfg["config_genome_size"] == "12"
        check cfg["config_block_size"] == "4096"

    #test "bad input":
    #    expect PbError:
    #        let cfg = parse("foo=bar")

suite "ipa_polish_prepare":
    test "countLines":
        let sample = """
0001
0002
"""
        let count = countLines(streams.newStringStream(sample))
        check count == 2
    test "getBlockIdFromFilename":
        check 0 == getBlockIdFromFilename("foo/bar.0.reads")
        check 99 == getBlockIdFromFilename("foobar.99.reads")

    test "shardUpperTriangular":
        proc sut(n, nshards: int): string =
            let shards = shardUpperTriangular(n = n, nshards = nshards)
            var outs = streams.newStringStream()
            for shard in shards:
                outs.writeLine(len(shard))
                for pr in shard:
                    outs.writeLine(pr)
            outs.setPosition(0)
            return outs.readAll()

        check sut(0, 0) == """
"""
        check sut(1, 1) == """
1
0 0 1
"""
        check sut(2, 1) == """
2
0 0 2
1 1 2
"""
        check sut(2, 2) == """
1
0 0 2
1
1 1 2
"""
        check sut(2, 3) == """
1
0 0 1
1
0 1 2
1
1 1 2
"""
        check sut(3, 1) == """
3
0 0 3
1 1 3
2 2 3
"""
        check sut(3, 2) == """
1
0 0 3
2
1 1 3
2 2 3
"""
        check sut(4, 1) == """
4
0 0 4
1 1 4
2 2 4
3 3 4
"""
        check sut(4, 2) == """
2
0 0 4
1 1 2
3
1 2 4
2 2 4
3 3 4
"""
        check sut(7, 1) == """
7
0 0 7
1 1 7
2 2 7
3 3 7
4 4 7
5 5 7
6 6 7
"""
        check sut(7, 2) == """
3
0 0 7
1 1 7
2 2 3
5
2 3 7
3 3 7
4 4 7
5 5 7
6 6 7
"""
        check sut(7, 3) == """
2
0 0 7
1 1 4
3
1 4 7
2 2 7
3 3 4
4
3 4 7
4 4 7
5 5 7
6 6 7
"""

    test "shardMatrix":
        proc shardm(nrows, ncols, nshards: int): string =
            let shards = shardMatrix(nrows = nrows, ncols = ncols, nshards = nshards)
            var outs = streams.newStringStream()
            for shard in shards:
                outs.writeLine(len(shard))
                for pr in shard:
                    outs.writeLine(pr)
            outs.setPosition(0)
            return outs.readAll()

        check shardm(0, 0, 0) == """
"""
        check shardm(1, 1, 1) == """
1
0 0 1
"""
        check shardm(2, 1, 1) == """
2
0 0 1
1 0 1
"""
        check shardm(2, 1, 2) == """
1
0 0 1
1
1 0 1
"""
        check shardm(1, 2, 1) == """
1
0 0 2
"""
        # With 3 shards, we snake back and forth.
        check shardm(2, 3, 3) == """
1
0 0 2
2
0 2 3
1 2 3
1
1 0 2
"""

    test "shardMatrixColumns":
        proc shardmc(nrows, ncols, nshards: int): string =
            let shards = shardMatrixColumns(nrows = nrows, ncols = ncols, nshards = nshards)
            var outs = streams.newStringStream()
            for shard in shards:
                outs.writeLine(len(shard))
                for pr in shard:
                    outs.writeLine(pr)
            outs.setPosition(0)
            return outs.readAll()

        check shardmc(0, 0, 0) == """
"""
        check shardmc(1, 1, 1) == """
1
0 0 1
"""
        check shardmc(2, 1, 1) == """
2
0 0 1
1 0 1
"""
        check shardmc(2, 1, 2) == """
2
0 0 1
1 0 1
"""
        check shardmc(1, 2, 1) == """
1
0 0 2
"""
        # With 3 shards, we snake back and forth.
        check shardmc(2, 3, 3) == """
2
0 0 1
1 0 1
2
0 1 2
1 1 2
2
0 2 3
1 2 3
"""

suite "ipa-separate-p-from-a":
    test "isHaplotigHeader":
        check isHaplotigHeader("foo") == false

suite "ipa2Header2Falcon":
    test "renamedSeq":
        expect system.Exception:
            discard renamedSeq("foo")
        check renamedSeq("foo.bar") == "bar"
        check renamedSeq("foo-bar") == "bar"
        check renamedSeq("foo.bar.baz") == ""
        check renamedSeq("foo-01-02") == ""
        check renamedSeq("foo-01-01") == "01_01"
