# vim: sts=4:ts=4:sw=4:et:tw=0

from streams import nil
from times import now, utc, `$`
import json
import unittest

# XrayOutputFormatter is very simple.
#
type
    XrayOutputFormatter* = ref object of unittest.OutputFormatter
        stream: streams.Stream
        startTime: times.DateTime
        #testErrors: seq[string]
        isInSuite: bool
        isInTest: bool
        failed: bool

proc newXrayOutputFormatter*(stream: streams.Stream): XrayOutputFormatter =
    return XrayOutputFormatter(
      stream: stream,
      startTime: now().utc,
      #testErrors: @[],
        failed: false,
    )

proc close*(formatter: XrayOutputFormatter) =
    ## Completes the report and closes the underlying stream.
    let endTime = now().utc
    let pf = if formatter.failed:
        "FAIL"
    else:
        "PASS"
    var xray_report = %* {
        "tests":
        [
            {
                "comment": "IPA falconc components",
                "finish": $endTime,
                "start": $formatter.startTime,
                "status": $pf,
                "testKey": "TAGT-448",
            }
        ]
    }
    streams.writeLine(formatter.stream, $xray_report)
    formatter.stream.close()

#method suiteStarted*(formatter: XrayOutputFormatter, suiteName: string) =
#  formatter.isInSuite = true

#method testStarted*(formatter: XrayOutputFormatter, testName: string) =
#  formatter.isInTest = true

method failureOccurred*(formatter: XrayOutputFormatter, checkpoints: seq[
        string], stackTrace: string) =
    #echo "FAILURE!"
    formatter.failed = true

method testEnded*(formatter: XrayOutputFormatter, testResult: TestResult) =
    formatter.isInTest = false
    var didNotFail = case testResult.status
        of OK: true
        of FAILED: false
        of SKIPPED: true
    if not didNotFail:
        #echo "FAIL!"
        formatter.failed = true

method suiteEnded*(formatter: XrayOutputFormatter) =
    formatter.isInSuite = false

######
