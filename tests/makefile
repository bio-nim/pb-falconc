MACPORTS_LIBDIR=/opt/local/lib
export DYLD_LIBRARY_PATH=${MACPORTS_LIBDIR}

all:
	rm -f ./all
	nim -d:nimDebugDlOpen c all.nim
	./all
each:
	nim c -r t_welcome.nim
	nim c -r t_kmers.nim
	nim c -r t_gc.nim
	nim c -r t_overlapFilter.nim
	nim c -r t_overlapParser.nim
	nim c -r t_ovl_cov_stats.nim
	nim c -r t_pbcromwell.nim
	nim c -r t_pbreports.nim
	nim c -r t_stats.nim
	nim c -r t_stats_gff.nim
	nim c -r t_rotate.nim
	nim c -r t_raptor_db.nim
	nim c -r t_align.nim
	nim c -r t_util.nim
	nim c -r t_ipa2.nim

.PHONY: all
