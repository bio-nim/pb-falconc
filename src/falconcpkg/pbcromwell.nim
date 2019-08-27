# vim: sw=4 ts=4 sts=4 tw=0 et:
from ./util import nil
from os import nil
from re import nil
from strutils import nil

let re_raw_long = re.re(r"\braw_reads\.\d+\.raw_reads\.\d+\.las$")
let re_p_long = re.re(r"\bpreads\.\d+\.preads\.\d+\.las$")

let re_raw_short = re.re(r"\braw_reads\.\d+\.las$")
let re_p_short = re.re(r"\bpreads\.\d+\.las$")

proc should_remove*(fn: string): bool =
    # For now, only check the base filename.
    if not strutils.endsWith(fn, ".las"):
        return false
    let basename = os.extractFilename(fn)
    if re.match(basename, re_raw_long):
        return true
    if re.match(basename, re_p_long):
        return true
    #let fi = os.getFileInfo(fn)
    return false

proc remove*(fn: string, dry_run: bool, verbose: bool = false) =
    # Assume file exists. (Otherwise might raise.)
    if verbose:
        if not dry_run:
            util.log("rm -f '" & fn & "'")
        else:
            util.log("#rm -f '" & fn & "'")
    if not dry_run:
        os.removeFile(fn) # can raise

proc remove_las*(run_dir: string, verbose: bool = false,
        dry_run: bool = false) =
    ## Remove all .las files except final stage of merge. Unzip is still possible.
    ## (Someday, we will add a flag to delete the final stage too, optionally.)
    if verbose:
        util.log("Deleting some las files under directory '" & run_dir & "' ...")
        if dry_run:
            util.log("(dry-run mode)")
    for fn in util.walk(run_dir, relative = true):
        if should_remove(fn):
            remove(fn, verbose = verbose, dry_run = dry_run)
