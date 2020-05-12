# vim: sw=4 ts=4 sts=4 tw=0 et:

import falconcpkg/ipa2_construct_config
import falconcpkg/ipa2_polish
import unittest, streams, tables

let expected_default_bash = """
config_autocomp_max_cov=1
config_block_size=4096
config_coverage=0
config_existing_db_prefix=''
config_genome_size=0
config_max_polish_block_mb=100
config_ovl_filter_opt='--max-diff 80 --max-cov 100 --min-cov 1 --bestn 10 --min-len 4000 --gapFilt --minDepth 4'
config_ovl_min_idt=98
config_ovl_min_len=1000
config_ovl_opt=''
config_phase_run=1
config_phasing_opt=''
config_phasing_piles=10000
config_phasing_split_opt='--split-type noverlaps --limit 3000000'
config_polish_run=1
config_seeddb_opt='-k 32 -w 80 --space 2'
config_seqdb_opt='--compression 0'
config_use_seq_ids=0
"""

let expected_default_json = """
{
    "config_autocomp_max_cov": "1",
    "config_block_size": "4096",
    "config_coverage": "0",
    "config_existing_db_prefix": "",
    "config_genome_size": "0",
    "config_max_polish_block_mb": "100",
    "config_ovl_filter_opt": "--max-diff 80 --max-cov 100 --min-cov 1 --bestn 10 --min-len 4000 --gapFilt --minDepth 4",
    "config_ovl_min_idt": "98",
    "config_ovl_min_len": "1000",
    "config_ovl_opt": "",
    "config_phase_run": "1",
    "config_phasing_opt": "",
    "config_phasing_piles": "10000",
    "config_phasing_split_opt": "--split-type noverlaps --limit 3000000",
    "config_polish_run": "1",
    "config_seeddb_opt": "-k 32 -w 80 --space 2",
    "config_seqdb_opt": "--compression 0",
    "config_use_seq_ids": "0"
}
"""

suite "ipa2_construct_config":
    test "write bash":
        var
            outs = streams.newStringStream()
            ins = streams.newStringStream()
        run(outs, ins, 'b', true)
        outs.setPosition(0)
        check outs.readAll() == expected_default_bash

    test "write json":
        var
            outs = streams.newStringStream()
            ins = streams.newStringStream()
        run(outs, ins, 'j', true)
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

    test "bad input":
        expect PbError:
            let cfg = parse("foo=bar")

suite "ipa_polish_prepare":
    test "countReadsInBlock":
        let sample = """
0001
0002
"""
        let count = countReadsInBlock(streams.newStringStream(sample))
        check count == 2
    test "getBlockIdFromFilename":
        check 0 == getBlockIdFromFilename("foo/bar.0.reads")
        check 99 == getBlockIdFromFilename("foobar.99.reads")
