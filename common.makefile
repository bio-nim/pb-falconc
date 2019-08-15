D=/pbi/flash/cdunn/bb/dipsim

default:
run:
	./falconc phasr -a $D/alignments/aln.0.001_0.001.sort.bam -r $D/rangen/random.fa -o results

pretty:
	find src integ-tests -name '*.nim' | xargs -L1 nimpretty --indent=4
