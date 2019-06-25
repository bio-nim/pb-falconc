# vim: sw=4 ts=4 sts=4 tw=0 et:
from os import nil
from pbpkg/welcome import nil

proc dataset(extras: seq[string]) =
    echo "pb dataset(" & $extras & ")"
    echo welcome.getWelcomeMessage()

when isMainModule:
    dataset(os.commandLineParams())
