SHELL=bash

D=/pbi/dept/secondary/siv/testdata/hgap/m4filt/greg200k-sv2

test:
	find $D -name '*.m4' | ../../src/falconc m4filt -n 1 --max-diff 120 --max-cov 120 --min-cov 2 --min-len 1 --out preads.ovl --filter-log filter.log
	diff preads.ovl $D/expect-preads.ovl
	diff filter.log $D/expect-filter.log
	rm -f preads.ovl filter.log

clean:
	rm -f preads.ovl filter.log
