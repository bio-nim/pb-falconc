I think this came from a run of greg200k-sv2.

To drop 98% of the lines:

    perl -lan -i -e 'print $_ if ((0 + "@F[1]") % 100 < 10 and (0 + "@F[0]") % 100 < 10)' preads.?.m4
