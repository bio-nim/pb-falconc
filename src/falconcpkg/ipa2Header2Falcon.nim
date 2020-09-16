from ./util import nil
from strformat import fmt
import hts
import strutils
import os

proc renameSeqs(seq_fn, output_prefix: string, pa: bool) =

    var refx: hts.Fai
    if not hts.open(refx, seq_fn):
        util.raiseEx(format("[FATAL] Could not open '$#'", seq_fn))

    let extension = (if pa : "a" else: "p")

    var f = open("{output_prefix}.{extension}.fasta".fmt, fmWrite)

    for i in 0 .. (refx.len - 1):
        let ctgSeq = refx.get(refx[i])
        let name_parts = refx[i].split({'.', '-'})
        var new_name: string
        if name_parts.len > 2:
            if name_parts[3] != "01":
                continue
            new_name = ">{name_parts[1]}_{name_parts[2]}".fmt
        else:
            new_name = ">{name_parts[1]}".fmt
        f.writeLine(new_name)
        f.writeLine(ctgSeq)
    f.close

proc main*(input_p_fn, input_a_fn, output_prefix: string) =
    ##Rename IPA2 fasta header names to match falcon
    if input_p_fn == "" or input_a_fn == "" or output_prefix == "":
        let msg = "Missing input or output required options."
        util.raiseEx(msg)
    renameSeqs(input_p_fn, output_prefix, false)
    renameSeqs(input_a_fn, output_prefix, true)

when isMainModule:
    main()
