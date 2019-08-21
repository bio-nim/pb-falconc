
import streams
import unittest

var f_of: File
if not open(f_of, "falconc_unit_tests_junit.xml", fmWrite):
    quit "Problem setting up report stream"
var fs_of = newFileStream(f_of)
var report_of = newJUnitOutputFormatter(fs_of)

#This allows the tests to be printed to the junit report
addOutputFormatter(report_of)

#This allows the tests to be printed to the console
addOutputFormatter(defaultConsoleFormatter())

#The includes need to happen after we setup the report formatters

include "t_gc.nim"
include "t_kmers.nim"
include "t_overlapFilter.nim"
include "t_pbcromwell.nim"
include "t_pbreports.nim"
include "t_stats.nim"
include "t_rotate.nim"

#cleanup
close(report_of)
close(fs_of)
