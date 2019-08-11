# vim: sw=4 ts=4 sts=4 tw=0 et:
#import falconcpkg/zev
from falconcpkg/align import nil
from falconcpkg/raptor_db import nil
from falconcpkg/rotate import nil
from falconcpkg/phasr import nil
from falconcpkg/overlapFilter import nil
from falconcpkg/pbcromwell import nil
from falconcpkg/pbreports import nil
from falconcpkg/stats import nil

const cToolVersion = "0.1.0"
const cGitSha1 {.strdefine.} = "0000000000000000000000000000000000000000"

proc version() =
    echo cToolVersion & "@" & cGitSha1
proc dataset(extras: seq[string]) =
    echo "falconc dataset"
proc kmers(int_dummy: int = 42, string_dummy: string = "hello") =
    echo "falconc kmers"
proc utils(extras: seq[string], float_req: float) =
    echo "falconc utils ..."
#proc zev() =
#    echo "starting"
#    zev.foo()
#    echo "finished"

when isMainModule:
    import cligen
    dispatchMulti(
        [version],
        [dataset, short = {}, help = {}],
        [kmers, short = {"int_dummy": 'd'}, help = {}],
        [utils, short = {}, help = {"float_req": "special help message"}],
        [align.align_filter, cmdName = "align-filter"],
        [raptor_db.filter, cmdName = "raptor-db-filter"],
        [rotate.main, cmdName = "circ-orient",
            help = {
       "input": "fasta file of circular sequences",
       "output": "fasta file output",
       "wl": "white list of sequences to rotate, one per line, no spaces, no trailing spaces",
       "window": "window size to caculate gc-skew",
       "step": "window step",
       "print": "print skew data to files, one per sequence"
        },
    ],
        [rotate.randomize, cmdName = "circ-randomize",
            help = {
       "input": "fasta file of circular sequences",
       "output": "fasta file output",
       "seed": "set seed, if non-zero",
        },
    ],
        [phasr.main, cmdName = "phasr",
          help = {
            "aln_fn": "BAM alignment, sorted on 'coordinate'",
            "ref_fn": "FASTA reference",
            "out_fn": "Output file name",
            "iterations": "Number of phasing iterations per read",
            "kmersize": "Kmer size",
            "delta": "Frequency cut",
        },
    ],
        [overlapFilter.runMergeBlacklists, cmdName = "m4filt-merge-blacklist",
        ],
        [overlapFilter.runStage1, cmdName = "m4filt-stage1",
        ],
        [overlapFilter.runStage2, cmdName = "m4filt-stage2",
        ],
        [overlapFilter.runDumpBlacklist, cmdName = "m4filt-dump-blacklist",
        ],
        [overlapFilter.falconRunner, cmdName = "m4filt-falconRunner",
         help = {
          "lasJson": "List of las files from falcon, e.g ../1-preads_ovl/las-merge-combine/las_fofn.json",
          "idtStage1": "Stage one percent identity filter, formatted as percentage, overlaps < %idt are skipped",
          "idtStage2": "Stage two percent identify filter",
          "minLen": "Minimum read length, reads shorter than minLen will be discarded",
          "minCov": "Minimum number of overlaps on either side of a read",
          "maxCov": "Maximum number of overlaps on either side of a read",
          "maxDiff": "Reads are skipped is abs(5p-3p) overlap counts > maxDiff",
          "bestN": "Keep N best overlaps at 5prime AND 3prime of a read",
          "minDepth": "Depths lower than minDepth are considered gaps",
          "gapFilt": "Run depth filter, takes a little more time",
          "nProc": "Number of processes to run locally",
          "filterLog": "Write read filter stats to this file",
          "outputFn": "Final m4 overlap file",
            }
        ],
        [overlapFilter.ipaRunner, cmdName = "m4filt-ipaRunner",
         help = {
          "ovlsFofn": "List of m4 files from ipa/raptor",
          "idtStage1": "Stage one percent identity filter, formatted as percentage, overlaps < %idt are skipped",
          "idtStage2": "Stage two percent identify filter",
          "minLen": "Minimum read length, reads shorter than minLen will be discarded",
          "minCov": "Minimum number of overlaps on either side of a read",
          "maxCov": "Maximum number of overlaps on either side of a read",
          "maxDiff": "Reads are skipped is abs(5p-3p) overlap counts > maxDiff",
          "bestN": "Keep N best overlaps at 5prime AND 3prime of a read",
          "minDepth": "Depths lower than minDepth are considered gaps",
          "gapFilt": "Run depth filter, takes a little more time",
          "nProc": "Number of processes to run locally",
          "filterLog": "Write read filter stats to this file",
          "outputFn": "Final m4 overlap file",
            }
        ],
        [pbcromwell.remove_las, cmdName = "pbcromwell-rm-las",
         short = {"dry_run" : 'n'},
        ],
        [pbreports.circ, cmdName = "pbreports-circ",
         help = {
          "fasta_fn": "FASTA filename, preferably indexed (with .fai)",
          "circ_fn": "Text list of circular contigs (newline-delimited)",
         }
        ],
        [stats.assembly, cmdName = "stats-assembly",
         help = {
          "fasta_fn": "FASTA filename, preferably indexed (with .fai)",
         }
        ],
    )
