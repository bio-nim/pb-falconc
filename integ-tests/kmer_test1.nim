# vim: sw=4 ts=4 sts=4 tw=0 et:
#from strformat import fmt
import deques
from os import nil
from strutils import format
from falconcpkg/kmers import Bin, difference, nkmers

proc log(label: string, formatstr: string, a: varargs[string, `$`]) =
    let msg = strutils.format(formatstr, a)
    #let tag = strutils.format("[TEST:$#] ", system.currentSourcePath(), label)
    let tag = strutils.format("[TEST:$#] ", label)
    echo tag, msg

proc main*(args: seq[string]): int =
    var sq = "ATCGGCTACTATT"

    echo format("Starting seq: $#", sq)

    var ans_lookups = [
        "AGCCGATGATAA",
        "TAGCCGATGATA",
        "ATCGGCTACTAT",
        "TCGGCTACTATT",
    ]
    var kms: kmers.pot_t = kmers.dna_to_kmers(sq, 12)
    var qms: kmers.pot_t = kmers.dna_to_kmers(sq, 12)

    echo "kms"
    kmers.print_pot(kms)
    echo ""
    echo "qms"
    kmers.print_pot(qms)

    var skms = kmers.initSpot(kms) # also sorts

    var final_res: int = 0
    var i: int32 = 0

    while i < skms.seeds.len():
        let tmp = kmers.bin_to_dna(skms.seeds[i].kmer, skms.word_size,
                                   skms.seeds[i].strand)
        let res = cmp(tmp, ans_lookups[i])

        log("DNA->BIT->DNA", "kmer: $# pos:$# expecting:$# observed:$# [$#]",
            skms.seeds[i].kmer, skms.seeds[i].pos, ans_lookups[i], tmp,
            if res != 0: "FAIL" else: "PASS")

        final_res = final_res or res

        inc(i)

    var hits = kmers.search(skms, qms)

    try:
        while true:
            let pair = deques.popLast(hits)
            echo format("qb:$# tb:$# qs:$# ts:$# $# $#",
                pair.a.pos, pair.b.pos,
                pair.a.strand, pair.b.strand,
                kmers.bin_to_dna(pair.a.kmer, skms.word_size, pair.a.strand),
                kmers.bin_to_dna(pair.b.kmer, skms.word_size, pair.b.strand))
    except:
        discard


    var tbins = [832088.Bin, 14315983.Bin, 3328355.Bin, 3578995.Bin]
    var fbins = [83.Bin, 0.Bin, 5.Bin, 10000.Bin]

    for b in tbins:
        let res = kmers.haskmer(skms, b)
        final_res = (not res).int
        log("HASKMER", "positive query:$# response:$# [$#]",
        b, res, if res != true: "FAIL" else: "PASS")

    for b in fbins:
        let res = kmers.haskmer(skms, b)
        final_res = final_res or res.int
        log("HASKMER", "negative query:$# response:$# [$#]",
         b, res, if res != false: "FAIL" else: "PASS")

    var sqms = kmers.initSpot(qms)
    var dqms = difference(skms, sqms)
    final_res = final_res or nkmers(dqms)

    return final_res

when isMainModule:
    quit main(os.commandLineParams())
