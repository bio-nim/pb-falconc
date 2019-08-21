# vim: sw=4 ts=4 sts=4 tw=0 et:
import falconcpkg/pbreports
import unittest
import json
import sets
import sequtils
import tables
from streams import nil

let expected_report = %*
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
                "values": [
                    "L2",
                    "C2",
                    "L1",
                    "C1"
                ]
            },
            {
                "header": "Length",
                "id": "microbial_asm_polishing_report.contigs_table.length",
                "values": [
                    25,
                    20,
                    15,
                    10
                ]
            },
            {
                "header": "Circular?",
                "id": "microbial_asm_polishing_report.contigs_table.circular",
                "values": [
                    false,
                    true,
                    false,
                    true
                ]
            }
        ],
        "id": "microbial_asm_polishing_report.contigs_table",
        "title": "Polished contigs from Microbial Assembly"
        }
    ],
    "tags": [],
    "title": "Microbial Assembly Polishing Report",
    "uuid": "???",
    "version": "1.0.0"
    }

suite "pbreports":
    test "get_report":
        let all_ctgs = {"C1": 10, "C2": 20, "L1": 15, "L2": 25}.toTable()
        let circ_ctgs = @["C1", "C2"]
        let j = pbreports.get_report(all_ctgs, circ_ctgs)
        j["uuid"] = % "???"
        check json.pretty(j) == json.pretty(expected_report)
