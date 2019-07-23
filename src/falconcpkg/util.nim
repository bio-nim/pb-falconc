# vim: sts=4:ts=4:sw=4:et:tw=0
import os

type PbError* = object of Exception

proc raiseEx*(msg: string) {.discardable.} =
    raise newException(PbError, msg)

proc isEmptyFile*(fin: string): bool =
    var finfo = getFileInfo(fin)
    if finfo.size == 0:
        return true
    return false

template withcd*(newdir: string, statements: untyped) =
    let olddir = os.getCurrentDir()
    os.setCurrentDir(newdir)
    defer: os.setCurrentDir(olddir)
    statements

proc log*(words: varargs[string, `$`]) =
    for word in words:
        write(stderr, word)
    write(stderr, '\l')
