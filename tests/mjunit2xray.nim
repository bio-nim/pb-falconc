
when isMainModule:
    ## This summarizes the multi-junit unit test into a single x-ray report.

    import parsexml
    import json
    import os
    import streams
    import strutils
    import strformat
    import times

    let starting_time = now().utc

    var failure = false

    if paramCount() < 1:
        quit("Usage: junit.xml")
    var x: XmlParser
    var fs = newFileStream(paramStr(1), fmRead)
    open(x, fs, paramStr(1))
    while true:
        x.next()
        case x.kind
        of xmlElementOpen:
            if cmpIgnoreCase(x.elementName, "failure") == 0:
                failure = true
        of xmlEof: break
        else: discard

    var pf = "PASS"
    if failure:
        pf = "FAIL"


    let ending_time = now().utc

    var xray_report = %* {"tests": [{"testKey": "TAGT-448",
            "start": $starting_time, "finish": $ending_time,
            "comment": "IPA falconc components", "status": "{pf}".fmt}]}
    var output = open("xray_summary.json", fmWrite)
    output.writeLine(xray_report)
    close(output)
