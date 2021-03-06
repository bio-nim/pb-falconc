# vim: sw=4 ts=4 sts=4 tw=0 et:
from kmers import hash  # avoiding "*" imports
import unittest
import deques
import sequtils
import sets

test "bin_to_dna":
    check kmers.bin_to_dna(0, 1, false) == "A"
    check kmers.bin_to_dna(1, 1, false) == "C"
    check kmers.bin_to_dna(2, 1, false) == "G"
    check kmers.bin_to_dna(3, 1, false) == "T"

    check kmers.bin_to_dna(0b00011011, 4, false) == "ACGT"
    check kmers.bin_to_dna(0b00011011, 4, true) == "TGCA"

test "dna_to_kmers":
    check kmers.dna_to_kmers("AAAA", 2).seeds.len() == 6

test "sorted_kmers":
    let
        sq = "ATCGGCTACTATT"
        expected = [
            "AGCCGATGATAA",
            "TAGCCGATGATA",
            "ATCGGCTACTAT",
            "TCGGCTACTATT",
        ]
        k = 12
    var
        kms = kmers.dna_to_kmers(sq, k)
    check kms != nil
    let spot = kmers.initSpot(kms) # sort
    check kms == nil
    let got = sequtils.mapIt(spot.seeds, kmers.bin_to_dna(it.kmer, k.uint8,
            it.strand))
    check got == expected
    #check kmers.haskmer("AGCCGATGATAA", kms)

test "search":
    let
        sq = "ATCGGCTACTATT"
        k = 12
        qms = kmers.dna_to_kmers(sq, k)
    var
        kms = kmers.dna_to_kmers(sq, k)
    let spot = kmers.initSpot(kms)
    let hits = kmers.search(spot, qms)
    check hits.len() == 4
    #check sets.toSet(seqUtils.toSeq(hits)).len() == 4 # 4 unique items
    check sets.len(sets.toSet(seqUtils.toSeq(deques.items(hits)))) == 4  # same as above

suite "difference":
    let
        sq = "ATCGGCTACTATT"
        k = 12

    test "difference_of_self_is_nothing":
        let kms = kmers.dna_to_kmers(sq, k)
        #let qms = deepCopy(kms)
        var qms: kmers.pot_t
        deepCopy(qms, kms)
        check qms[] == kms[]
        check kmers.nkmers(qms) == 4
        check kmers.nkmers(kms) == 4
        let qspot = kmers.initSpot(qms)
        check qms == nil

        let kms0 = kmers.difference(kms, qspot)

        check kmers.nkmers(kms0) == 0
        check kmers.nkmers(kms) == 4
        check kmers.nkmers(qspot) == 4

        let
            expected: array[0, string] = []
            got = kmers.get_dnas(kms0)
        check got == expected

    test "difference_of_nothing_is_self":
        let kms = kmers.dna_to_kmers(sq, k)
        var qms = kmers.dna_to_kmers("", k)
        #let orig = deepCopy(kms)
        var orig: kmers.pot_t
        deepCopy(orig, kms)
        check kmers.nkmers(qms) == 0
        check kmers.nkmers(kms) == 4
        let qspot = kmers.initSpot(qms)
        check qms == nil

        let kms4 = kmers.difference(kms, qspot)

        check kmers.nkmers(kms4) == 4
        check kmers.nkmers(kms) == 4
        let got = kmers.get_dnas(kms4)
        let expected = kmers.get_dnas(orig)
        check got == expected
