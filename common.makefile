D=/pbi/flash/cdunn/bb/dipsim

default:
run:
	./falconc phasr -a $D/alignments/aln.0.001_0.001.sort.bam -r $D/rangen/random.fa -o results

pretty:
	#find src tests integ-tests -name '*.nim' | xargs -L1 /pbi/flash/cdunn/gh/Nim/bin/nimpretty --indent=4 --maxLineLen:200
	find src tests integ-tests -name '*.nim' | xargs -L1 nimpretty --indent=4 --maxLineLen:200
