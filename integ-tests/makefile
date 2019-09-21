default: test1 test-phasr
test1:
	nim c kmer_test1
	./kmer_test1 >| kmer_test1.result.txt
	diff kmer_test1.result.txt kmer_test1.expect.txt
test-phasr:
	nim c phasr_test
	./phasr_test
alarm:
	# run from top-dir
	-src/falconc raptor-db-calc-length-cutoff -r ../raptor/test-data/raptordb-fetch/test-2.rdb --genome 10000000 -f --alarms foo
	cat foo
	rm -f foo
