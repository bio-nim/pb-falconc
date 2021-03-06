# vim: sw=4 ts=4 sts=4 tw=0 et:
import strutils
import hts

from ./align import calc_query_pos
from ./util import raiseEx
from ./gc import revComp


proc bam_to_fasta*(in_bam, region: string, flag: int = 3844,
        flip_rc: bool = false) =
    ## Very similar to `samtools fasta`, but the reads are subSequenced.  Output is printed to STDOUT

    let region = region.split({':'});

    if region.len != 2:
        raiseEx(format("[FATAL] Region is incorrectly formatted '$#'", region.join(" ")))

    let pos    = region[1].split({'-'})

    if pos.len != 2:
        raiseEx(format("[FATAL] Region is incorrectly formatted '$#'", pos.join(" ")))

    var b: Bam
    open(b, in_bam, index = true)
    defer:
        b.close()



    let t_start = parseInt(pos[0])
    let t_end = parseInt(pos[1])


    for r in b.query(region[0], t_start, t_end):
        if ((r.flag.int and flag) > 0):
            continue

        let qpos = align.target_to_query(r.cigar, r.start, t_start, t_end)

        var dna_seq: string
        r.sequence(dna_seq)
        dna_seq = dna_seq[qpos.q_start .. qpos.q_end]
        echo ">", r.qname, " len:", dna_seq.len
        if flip_rc and r.flag.reverse:
            dna_seq = gc.revComp(dna_seq)
        echo dna_seq
