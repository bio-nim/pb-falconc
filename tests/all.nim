
import streams
import unittest

var fs_of = streams.newFileStream("falconc_unit_tests_junit.xml", fmWrite)
var report_of = unittest.newJUnitOutputFormatter(fs_of)

#This allows the tests to be printed to the junit report
unittest.addOutputFormatter(report_of)

#This allows the tests to be printed to the console
unittest.addOutputFormatter(unittest.defaultConsoleFormatter())

#The includes need to happen after we setup the report formatters

include "t_gc.nim"
include "t_kmers.nim"
include "t_overlapFilter.nim"
include "t_pbcromwell.nim"
include "t_pbreports.nim"
include "t_stats.nim"
include "t_rotate.nim"

#cleanup
unittest.close(report_of)
streams.close(fs_of)
