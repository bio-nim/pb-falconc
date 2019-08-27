# vim: sw=4 ts=4 sts=4 tw=0 et:
from ./util import nil
from ./nuuid import nil
from algorithm import nil
from sets import nil
from streams import nil
from tables import `$`, `[]=`, pairs
from hts import `[]`
import json

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

proc get_report*(unsorted_columns: seq[Column]): json.JsonNode =
    # Sort in descending order of length.
    let columns = algorithm.sortedByIt(unsorted_columns, -it.length)

    # Unzip.
    let n = columns.len()
    var
        values_contig = newSeqOfCap[string](n)
        values_length = newSeqOfCap[int](n)
        values_circular = newSeqOfCap[bool](n)
    for col in columns:
        values_contig.add(col.contig)
        values_length.add(col.length)
        values_circular.add(col.circular)

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
                        "header": "Circular?",
                        "id": "microbial_asm_polishing_report.contigs_table.circular",
                        "values": values_circular
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

proc get_report*(all_ctgs: tables.Table[string, int], circ_ctgs: seq[
        string]): json.JsonNode =
    var columns = newSeq[Column]()
    let circs = sets.toHashSet[string](circ_ctgs)
    for ctg, length in tables.pairs(all_ctgs):
        let isCirc = sets.contains(circs, ctg)
        let column: Column = (contig: ctg, length: length, circular: isCirc)
        columns.add(column)
    return get_report(columns)

proc circ*(fasta_fn: string, circ_fn: string) =
    ## Given FASTA of all ctgs and text-list of circular ctgs,
    ## print a report (pbreports format).
    let all_ctgs = get_all_ctgs(fasta_fn)
    let circ_ctgs = get_circ_ctgs(circ_fn)
    let j = get_report(all_ctgs, circ_ctgs)
    echo json.pretty(j, indent = 4)
