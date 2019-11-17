# vim: sw=4 ts=4 sts=4 tw=0 et:
import sets
import strutils
import hts
import strformat

proc parsePhaseInfo(phaseInfo: string): sets.HashSet[string] =
    ## parses the reads2ctg file, skipping the first line (comment).

    let f = open(phaseInfo, fmRead)
    defer: f.close()

    for line in f.lines:
        var lineDat = line.split()
        if lineDat[1] == "ctg":
            continue
        # unphased reads should be skipped - we will look them up later.
        if parseInt(lineDat[2]) == -1:
            continue
        sets.incl(result, lineDat[0])


proc getUnPhasedAssignment(bam_fn: string, skip: sets.HashSet[string]): seq[
  string] =
    ## Reads the bam file of unphased reads and pulls their assignment by primary mapping.
    var b: Bam
    hts.open(b, bam_fn, index = false)
    defer: b.close()

    let h = hts.targets(b.hdr)

    result = newSeq[string]()

    var n = 0
    var s = 0

    for record in b:
        inc(n)
        if record.qname in skip:
            inc(s)
            continue
        # Terrible mapping, skip. HARDCODED
        if record.mapping_quality < 1:
            inc(s)
            continue
        # Check sam flags, skip. HARDCODED
        if (record.flag.int and 0b11100000100) > 0: # 1796
            inc(s)
            continue
        result.add("{record.qname} {h[record.tid].name} -1 -1 -1".fmt)

proc runner*(phase_fn, bam_fn, out_fn: string) =
    ## Assigns the unphased reads to haplotigs or primary contigs based on
    ## pbmm2 mapping.

    let phasedReads = parsePhaseInfo(phase_fn)
    let toAdd = getUnPhasedAssignment(bam_fn, phasedReads)

    let fin = open(phase_fn, fmRead)
    defer: fin.close()

    let fout = open(out_fn, fmWrite)
    defer: fout.close()

    for line in fin.lines:
        var lineDat = line.split()
        if lineDat[1] == "ctg":
            fout.writeLine(line);
            continue
        if parseInt(lineDat[2]) == -1:
            continue
        fout.writeLine(line);

    for augs in toAdd:
        fout.writeLine(augs)
