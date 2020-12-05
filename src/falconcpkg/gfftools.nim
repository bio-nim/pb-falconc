# vim: sw=4 ts=4 sts=4 tw=0 et:
import tables
from algorithm import nil
from parsecsv import nil
from sequtils import nil
from streams import nil
from strformat import fmt
from strutils import startsWith
from ./util import PbError

type
    GffRow* = object
        sequence: string
        source: string
        feature: string
        tstart: int
        tend: int
        score: string
        strand: char
        phase: char
        attributes: string
    GffLine* = array[9, string]
    GffLines* = Table[string, seq[GffRow]]

var
    SEP* = '\t'

proc `$`(row: GffRow): string =
    return "{row.sequence}{SEP}{row.source}{SEP}{row.feature}{SEP}{row.tstart}{SEP}{row.tend}{SEP}{row.score}{SEP}{row.strand}{SEP}{row.phase}{SEP}{row.attributes}".fmt

proc loadGffLines*(sin: streams.Stream, fn: string, sep: char): ref GffLines =
    # Allow 8-field line instead of 9.
    new(result)
    var
        csv: parsecsv.CsvParser
        row: GffRow
    parsecsv.open(csv, sin, fn, separator = sep)

    while parsecsv.readRow(csv):
        if csv.row.len() == 0:
            continue
        elif csv.row[0].startsWith("#"):
            continue
        elif csv.row.len() != 9:
            raise newException(PbError, $csv.row)
        row.sequence = csv.row[0]
        row.source = csv.row[1]
        row.feature = csv.row[2]
        row.tstart = strutils.parseInt(csv.row[3])
        row.tend = strutils.parseInt(csv.row[4])
        row.score = csv.row[5]
        row.strand = csv.row[6][0]
        row.phase = csv.row[7][0]
        row.attributes = csv.row[8]
        let name = row.sequence
        if not result.contains(name):
            result[name] = @[row]
        else:
            result[name].add(row)

proc loadGffLines*(sin: streams.Stream, fn = "<stream>"): ref GffLines =
    # Allow space delimiters instead of tab.
    let start = streams.getPosition(sin)
    try:
        return loadGffLines(sin, fn, ' ')
    except PbError as e:
        streams.setPosition(sin, start)  # rewind
        return loadGffLines(sin, fn, '\t')

proc gffsubtractStreams*(gsin, masksin, sout: streams.Stream) =
    var
        glines = loadGffLines(gsin)
        masklines = loadGffLines(masksin)

    proc regionCmp(a,b: GffRow): int =
        # in perl, ($a[3] <=> $b[3] or $a[4] <=> $b[4])
        let cmpstart = system.cmp(a.tstart, b.tstart)
        if cmpstart != 0:
            return cmpstart
        let cmpend = system.cmp(a.tend, b.tend)
        return cmpend

    for seqname in algorithm.sorted(sequtils.toSeq(tables.keys(glines))):
        if seqname notin masklines:
            continue
        var
            gff2s = algorithm.sorted(masklines[seqname], regionCmp)
        for gff1 in glines[seqname]:
            var
                start = gff1.tstart
            for gff2 in gff2s:
                if gff2.tend < gff1.tstart:
                    continue
                if gff2.tstart > gff1.tend:
                    break
                if gff2.tstart > start:
                    var newrow = gff1
                    newrow.tstart = start
                    newrow.tend = gff2.tstart - 1
                    streams.writeLine(sout, $newrow)
                start = gff2.tend + 1
            if start <= gff1.tend:
                var newrow = gff1
                newrow.tstart = start
                streams.writeLine(sout, $newrow)

proc gffsubtract*(gff_fn, mask_gff_fn: string) =
    ## cut out everything in one GFF file that overlaps with a second

    # From:
    #   /mnt/software/g/gfftools/dalexander/gffsubtract.pl
    # which prints slightly more than the original:
    #   https://github.com/ihh/gfftools/blob/master/gffsubtract.pl

    # Print the headers of the 1st one
    for line in lines(gff_fn):
        if line.startsWith("##"):
            echo line
    echo "##mask {mask_gff_fn}\n".fmt
    var
        gsin = streams.openFileStream(gff_fn)
        msin = streams.openFileStream(mask_gff_fn)
    gffsubtractStreams(gsin, msin, streams.newFileStream(stdout))
