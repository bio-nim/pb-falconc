# vim: sw=4 ts=4 sts=4 tw=0 et:
from ./util import isEmptyFile, log
from os import quoteShell

proc run*() =
    echo "run"

proc main*(out_fn: string, out_fmt: string = "json") =
    ## Takes an advanced options string, and reformats it into JSON format.
    ## Input/output is on stdin/stdout. Options which aren't set explicitly in the input
    ## will be set to default (configurable via args).

    echo "Hi runner", out_fn, out_fmt
