# vim: sw=4 ts=4 sts=4 tw=0 et:
import os

proc dataset(extras: seq[string]) =
    echo "pb dataset(" & $extras & ")"

when isMainModule:
    dataset(os.commandLineParams())
