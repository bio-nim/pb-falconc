# vim: sw=4 ts=4 sts=4 tw=0 et:
from os import nil
from falconcpkg/welcome import nil

proc dataset(extras: seq[string]) =
    echo "falconc dataset(" & $extras & ")"
    echo welcome.getWelcomeMessage()

when isMainModule:
    dataset(os.commandLineParams())
