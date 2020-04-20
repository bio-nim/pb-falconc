# vim: sw=4 ts=4 sts=4 tw=0 et:
from ./util import isEmptyFile, log, raiseEx, PbError
from os import quoteShell
from strformat import `&`
import json, tables, sets, sequtils, strutils, streams

export PbError

var default_config = """
config_genome_size = 0
config_coverage = 0
config_polish_run = 1
config_phase_run = 1

config_existing_db_prefix =
config_block_size = 4096
config_seqdb_opt = --compression 0
config_seeddb_opt = -k 32 -w 80 --space 2
config_ovl_opt =
config_ovl_min_idt = 98
config_ovl_min_len = 1000
config_autocomp_max_cov = 1
config_ovl_filter_opt = --max-diff 80 --max-cov 100 --min-cov 1 --bestn 10 --min-len 4000 --gapFilt --minDepth 4
config_use_seq_ids = 0
config_phasing_opt =
config_phasing_piles = 10000
config_max_polish_block_mb = 100
"""

type
    ConfigTable* = tables.OrderedTable[string, string]
    #ConfigTuple = tuple
    #  genome_size,
    #    coverage,
    #    polish_run,
    #    phase_run,
    #    existing_db_prefix,
    #    block_size,
    #    seqdb_opt,
    #    seeddb_opt,
    #    ovl_opt,
    #    ovl_min_idt,
    #    ovl_min_len,
    #    ovl_filter_opt,
    #    use_seq_ids,
    #    phase_opt,
    #    phasing_piles: string

proc formatter_json(config_dict: ConfigTable): string =
    return json.pretty( % config_dict, 4)

proc formatter_bash(config_dict: ConfigTable): string =
    var lines: seq[string]
    for key in tables.keys(config_dict):
        let val = quoteShell(config_dict[key])
        # Accumulate the lines.
        let line = &"{key}={val}"
        lines.add(line)
    let ret = strutils.join(lines, "\n")
    return ret

proc convert_config_to_dict*(config: string): ConfigTable =
    # Replace newlines with a ';' to make it uniform. (Also to support ';' as separator.)
    let in_config = strutils.replace(config, '\n', ';')

    for line in strutils.split(in_config, ';'):
        let stripped = strutils.strip(line)
        if stripped == "": continue
        let keyval = strutils.split(stripped, '=')
        doAssert len(keyval) == 2, "Malformed config option. Each config option needs to have exactly one '=' character, specifying the config option on the left side, and the value on the right side. Line: " & $stripped
        let
            param_name = strutils.strip(keyval[0])
            param_val = strutils.strip(keyval[1])
        result[param_name] = param_val

proc validate_param_names(valid_params, test_params: ConfigTable) =
    let valid = sets.toHashSet(sequtils.toSeq(tables.keys(valid_params)))
    var unknown_params: seq[string]
    for p in tables.keys(test_params):
        if not sets.contains(valid, p):
            unknown_params.add(p)
    if len(unknown_params) != 0:
        let msg = "Unknown config parameters specified: " & $unknown_params
        raiseEx(msg)

proc parse*(in_str: string): ConfigTable =
    # Load the defaults.
    let default_config_dict = convert_config_to_dict(default_config)

    # Load the user specified options.
    let user_config_dict = convert_config_to_dict(in_str)

    # Validate the user config. There shouldn't be any keys which do not
    # appear in the default_config_dict.
    validate_param_names(default_config_dict, user_config_dict)

    # Update the defaults with user specified values.
    result = default_config_dict
    for k, v in tables.pairs(user_config_dict):
        result[k] = v

proc run*(fp_out, fp_in: streams.Stream, out_fmt: char, sort: bool) =
    # Collect the input lines.
    let in_str = strutils.join(sequtils.toSeq(lines(fp_in)), ";")

    var config_dict = parse(in_str)
    if sort:
        # Sort to match Python behavior.
        tables.sort(config_dict, system.cmp)

    # Write the dict.
    let formatter = if out_fmt == 'j': formatter_json else: formatter_bash
    let out_str = formatter(config_dict)
    fp_out.writeLine(out_str)

proc main*(out_fn: string, out_fmt: string = "json", in_fn = "-", no_sort = false) =
    ## Takes an advanced options string, and reformats it into JSON format.
    ## Input/output is on stdin/stdout. Options which aren't set explicitly in the input
    ## will be set to default (configurable via args).

    var fp_in = if in_fn == "-": streams.newFileStream(stdin) else: streams.newFileStream(in_fn, fmRead)
    var fp_out = streams.newFileStream(out_fn, fmWrite)
    run(fp_out, fp_in, out_fmt[0], not no_sort)
    fp_out.close()
