from ./util import nil
from strformat import fmt
import hts
import strutils

proc renamedSeq*(name: string): string =
    let name_parts = name.split({'.', '-'})
    if name_parts.len > 2:
        if name_parts[2] != "01":
            return ""
        return "{name_parts[1]}_{name_parts[2]}".fmt
    else:
        return "{name_parts[1]}".fmt

proc renameSeqs(seq_fn, output_prefix: string, pa: bool) =
    var refx: hts.Fai
    if not hts.open(refx, seq_fn):
        util.raiseEx(format("[FATAL] Could not open '$#'", seq_fn))

    let extension = (if pa : "a" else: "p")

    var f = open("{output_prefix}.{extension}.fasta".fmt, fmWrite)

    for i in 0 .. (refx.len - 1):
        let ctgSeq = refx.get(refx[i])
        let new_name = renamedSeq(refx[i])
        if new_name == "":
            continue
        f.write('>')
        f.writeLine(new_name)
        f.writeLine(ctgSeq)
    f.close

proc main*(input_p_fn, input_a_fn, output_prefix: string) =
    ##Rename IPA2 fasta header names to match falcon
    renameSeqs(input_p_fn, output_prefix, false)
    renameSeqs(input_a_fn, output_prefix, true)
