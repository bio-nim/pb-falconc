# vim: sw=4 ts=4 sts=4 tw=0 et:
from ./util import nil
from ./nuuid import nil
from algorithm import nil
from sets import nil
from streams import nil
from hts import `[]`
import tables
import json
import strutils
import stats_gff

proc get_all_ctgs(fasta: hts.Fai): tables.Table[string, int] =
    ## Given indexed FASTA, return all name->length.
    # (Hard to test, but very fast, since htslib can fseek to each header.)
    let n = hts.len(fasta)
    result = tables.initTable[string, int]()
    for i in 0 ..< n:
        let name = fasta[i]
        result[name] = hts.chrom_len(fasta, name)

proc get_all_ctgs*(fasta_fn: string): tables.Table[string, int] =
    var idx: hts.Fai
    if not hts.open(idx, fasta_fn):
        util.raiseEx("Could not open FASTA '" & fasta_fn & "'")
    return get_all_ctgs(idx)

proc get_circ_ctgs*(sin: streams.Stream): seq[string] =
    ## Simple textfile reader, one ctg per line
    var ctg: string
    while streams.readLine(sin, ctg):
        add(result, ctg)

proc get_circ_ctgs(collected_circ_txt_fn: string): seq[string] =
    let snarfed = system.readFile(collected_circ_txt_fn) # fast
    var sin = streams.newStringStream(snarfed)
    defer: streams.close(sin)
    return get_circ_ctgs(sin)

type
    Column = tuple
        contig: string
        length: int
        circular: bool
        coverage: float64

proc get_contig_table_report*(unsorted_columns: seq[Column]): json.JsonNode =
    # Sort in descending order of length.
    let columns = algorithm.sortedByIt(unsorted_columns, -it.length)

    # Unzip.
    let n = columns.len()
    var
        values_contig = newSeqOfCap[string](n)
        values_length = newSeqOfCap[int](n)
        values_circular = newSeqOfCap[bool](n)
        values_coverage = newSeqOfCap[float64](n)
    for col in columns:
        values_contig.add(col.contig)
        values_length.add(col.length)
        values_circular.add(col.circular)
        values_coverage.add(col.coverage)

    let uuid = nuuid.generateUUID()
    let version = "1.0.0"
    result = %*
        {
            "_comment": "pbreports-style JSON",
            "attributes": [],
            "dataset_uuids": [],
            "id": "microbial_asm_polishing_report",
            "plotGroups": [],
            "tables": [
                {
                "columns": [
                    {
                        "header": "Contig",
                        "id": "microbial_asm_polishing_report.contigs_table.contig",
                        "values": values_contig
                    },
                    {
                        "header": "Length",
                        "id": "microbial_asm_polishing_report.contigs_table.length",
                        "values": values_length
                    },
                    {
                        "header": "Circular",
                        "id": "microbial_asm_polishing_report.contigs_table.circular",
                        "values": values_circular
                    },
                    {
                        "header": "Coverage",
                        "id": "microbial_asm_polishing_report.contigs_table.coverage",
                        "values": values_coverage
                    }
                ],
                "id": "microbial_asm_polishing_report.contigs_table",
                "title": "Polished contigs from Microbial Assembly"
                },
            ],
            "tags": [],
            "title": "Microbial Assembly Polishing Report",
            "uuid": uuid,
            "version": version
        }

proc get_contig_table_report*(all_ctgs: tables.Table[string, int], circ_ctgs: seq[
        string], cov_dict: Table[string, Table[string, string]]): json.JsonNode =
    var columns = newSeq[Column]()
    let circs = sets.toHashSet[string](circ_ctgs)
    for ctg, length in tables.pairs(all_ctgs):
        let isCirc = sets.contains(circs, ctg)
        # Fetch the coverage. Alignment is to the draft, so we need to get
        # rid of the Arrow suffix from the contig name first.
        let ctgBasename = split(ctg, '|')
        # Allow nonexistent keys. Small and bad contigs could have zero entries in the GFF.
        var ctgCov = 0.0
        if ctg_basename[0] in cov_dict:
            ctgCov = parseFloat(cov_dict[ctg_basename[0]]["cov_mean"])
        # Form a column.
        let column: Column = (contig: ctg, length: length, circular: isCirc, coverage: ctgCov)
        columns.add(column)
    return get_contig_table_report(columns)

proc circ*(fasta_fn: string, circ_fn: string) =
    ## Given FASTA of all ctgs and text-list of circular ctgs,
    ## print a report (pbreports format).
    let all_ctgs = get_all_ctgs(fasta_fn)
    let circ_ctgs = get_circ_ctgs(circ_fn)
    # For legacy purposes, generate a dummy coverage dict.
    let cov_dict = initTable[string, Table[string, string]]()
    let j = get_contig_table_report(all_ctgs, circ_ctgs, cov_dict)
    echo json.pretty(j, indent = 4)

proc run_gen_contig_table*(fasta_fn: string, circ_fn: string, gff_fn: string) =
    # Parse contigs and lengths.
    let all_ctgs = get_all_ctgs(fasta_fn)

    # Parse the list of circular contigs.
    let circ_ctgs = get_circ_ctgs(circ_fn)

    # Parse the GFF file.
    let stream_gff = streams.newFileStream(gff_fn, fmRead)
    defer: stream_gff.close()
    let cov_dict = summarize_gff_coverage(stream_gff)

    # Generate the report.
    let j = get_contig_table_report(all_ctgs, circ_ctgs, cov_dict)
    echo json.pretty(j, indent = 4)

