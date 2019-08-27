# vim: sw=4 ts=4 sts=4 tw=0 et:
from ./util import nil
from os import nil
from re import nil
from strformat import nil
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
    return false

proc remove*(fn: string, dry_run: bool, verbose: int): BiggestInt =
    # Return size that is (or would be) removed.
    # Assume file exists. (Otherwise might raise.)
    if verbose > 1:
        if not dry_run:
            util.log("rm -f '" & fn & "'")
        else:
            util.log("#rm -f '" & fn & "'")
    #let fi = os.getFileInfo(fn, followSymlink=false)
    #result = fi.size
    if not dry_run:
        os.removeFile(fn) # can raise

proc remove_las*(run_dir: string, verbose: int = 1,
        dry_run: bool = false) =
    ## Remove all .las files except final stage of merge. Unzip is still possible.
    ## (Someday, we will add a flag to delete the final stage too, optionally.)
    if verbose > 0:
        util.log("Deleting some las files under directory '" & run_dir & "' ...")
        if dry_run:
            util.log("(dry-run mode)")
    var total: BiggestInt = 0
    var nFound = 0
    var nRemoved = 0
    var reportAt = 1
    util.withcd(run_dir):
        #for fn in util.walk(".", relative = true):
        for fn in lines(stdin):
            echo fn
            if true:
                continue
            nFound += 1
            if should_remove(fn):
                total += remove(fn, verbose = verbose, dry_run = dry_run)
                nRemoved += 1
            if nFound == reportAt:
                if verbose > 0:
                    reportAt *= 2
                    let msg = strformat.fmt(
                        "Found={nFound:8d}, removed={nRemoved:8d}, bytes={total:12d}, current: {fn}")
                    util.logt(msg)
    if verbose > 0:
        let msg = strformat.fmt(
            "Found={nFound:8d}, removed={nRemoved:8d}, bytes={total:12d}")
        util.logt(msg)
