# vim: sw=4 ts=4 sts=4 tw=0 et:
from strformat import fmt

proc split*(max_nshards: int, read_to_contig_fn, out_shard_prefix, out_block_prefix: string) =
    echo "prepare {read_to_contig_fn} {max_nshards} {out_shard_prefix} {out_block_prefix}".fmt
