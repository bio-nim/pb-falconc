# vim: sts=2:ts=2:sw=2:et:tw=0
from streams import nil
import unittest
import ./xray

var fs_of = streams.newFileStream("falconc_unit_tests_junit.xml", fmWrite)
var report_of = unittest.newJUnitOutputFormatter(fs_of)

#This allows the tests to be printed to the junit report
unittest.addOutputFormatter(report_of)

#This allows the tests to be printed to the console
unittest.addOutputFormatter(unittest.defaultConsoleFormatter())

#For importing into X-ray
var xfs_of = streams.newFileStream("xray_summary.json", fmWrite)
var xray_of = newXrayOutputFormatter(xfs_of)
unittest.addOutputFormatter(xray_of)


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
close(xray_of)
