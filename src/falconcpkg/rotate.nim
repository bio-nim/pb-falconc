import hts
import gc
import strutils
import algorithm
import random
import sets
import std/wordwrap
import logging
import os


var logger* = newConsoleLogger(fmtStr = "[$time] - $levelname: ")
addHandler(logger)

type

    dat* = tuple[pos: int, skew: float, accum: float]

    skewDat* = object
        data*: seq[dat]
        mini*: int
        minv*: float

proc checkEmptyFile(fin: string): bool =
    var finfo = getFileInfo(fin)
    if finfo.size == 0:
        return true
    return false

proc calcSkew*(dna: string, win: int, step: int): skewDat =
    var sk: skewDat
    var accum: float = 0
    var i, j: int = 0;
    var min: float = 100000000

    while(i <= dna.len - win - 1):
        var c = gc.countBases(dna[i..i+win])
        var d: dat
        var skew = gc.gcSkew(c)

        accum += skew
        d.pos = i
        d.skew = skew
        d.accum = accum
        sk.data.add(d)

        if accum < min:
            sk.mini = j
            sk.minv = accum
            min = accum

        inc(j)
        i += step

    return sk

proc printSkew*(seqname: string, sk: skewDat) =
    var output = open(seqname & ".gc_skew.txt", fmWrite)
    for i in sk.data:
        output.write(seqname, " ", i.pos, " ", i.skew.formatFloat(ffDecimal,
                2), " ", i.accum.formatFloat(ffDecimal, 2), "\n")
    output.close()

proc randomize*(input: string, output: string, seed: int64 = 0) =
    ##randomly rotates left rotates the sequence and writes to the output file.
    if seed != 0:
        random.randomize(seed)
    else:
        random.randomize()
    var fai: Fai
    if not fai.open(input):
        quit "couldn't open fasta"

    var outfh = open(output, fmWrite)
    defer:
        close(outfh)
    for i in 0..<fai.len:
        let chrom_name = fai[i]
        let chrom_len = fai.chrom_len(chrom_name)
        var full_sequence = fai.get(chrom_name)
        var ran = rand(chrom_len)
        discard algorithm.rotateLeft(full_sequence, ran)
        echo "[INFO] randomizing: ", chrom_name
        outfh.write(">", chrom_name, " randomly_shifted_by_bp:-",
                 ran, "/", chrom_len, "\n")
        outfh.write(wrapWords(full_sequence), "\n")


proc loadWhiteList*(fin: string): HashSet[string] =
    var whitelist = initHashSet[string]()
    if fin == "":
        return whitelist
    for l in fin.lines:
        echo l
        whitelist.incl(l)
    return whitelist

proc reorient(fin: string, fon: string, wl: string, w: int, s: int,
        print: bool) =

    var output = open(fon, fmWrite)
    defer:
        output.close()

    if checkEmptyFile(fin):
        logger.log(lvlNotice, "Empty input, output will be empty.")
        return

    var fai: Fai
    if not fai.open(fin):
        logger.log(lvlFatal, "Problem loading fasta file")
        quit 1

    var whiteList = loadWhiteList(wl)

    for i in 0..<fai.len:
        var chrom_name = fai[i]
        var chrom_len = fai.chrom_len(chrom_name)
        var full_sequence = fai.get(chrom_name)
        var sdf = calcSkew(full_sequence, w, s)
        if not (whiteList.contains(chrom_name)) and (whiteList.len > 0):
            sdf.data[sdf.mini].pos = 0
        if print:
            printSkew(chrom_name, sdf)
        echo "#windows:", sdf.data.len, " pivot index:", sdf.mini
        echo "pivot pos: ", sdf.data[sdf.mini].pos, " / ", chrom_len,
         " seq: ", chrom_name
        discard algorithm.rotateLeft(full_sequence, sdf.data[sdf.mini].pos)
        output.write(">", chrom_name, " shifted_by_bp:-",
                   sdf.data[sdf.mini].pos, "/", chrom_len, "\n")
        output.write(wrapWords(full_sequence), "\n")


proc main*(input_fn: string, output_fn: string, wl_fn = "", window = 500, step = 200,
        print = false) =
    ##reorients circular sequences based on gc-skew distribution and writes to output.
    if input_fn == "" or output_fn == "":
        logger.log(lvlFatal, "Missing input or output required options.")
        quit 1
    logger.log(lvlInfo, "Reorienting.")
    reorient(input_fn, output_fn, wl_fn, window, step, print)

when isMainModule:
    main()
