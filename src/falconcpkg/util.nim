# vim: sts=4:ts=4:sw=4:et:tw=0
#from cpuinfo import nil
from os import nil
#from threadpool import nil
from streams import nil
from strformat import nil
import osproc
import times

type PbError* = object of Exception
type GenomeCoverageError* = object of PbError
type FieldTooLongError* = object of PbError
type TooFewFieldsError* = object of PbError

proc raiseEx*(msg: string) {.discardable.} =
    raise newException(PbError, msg)

proc isEmptyFile*(fn: string): bool =
    var finfo = os.getFileInfo(fn)
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

proc logt*(words: varargs[string, `$`]) =
    var then {.global.} = times.now()
    let
        since = times.initDuration(seconds = times.inSeconds(times.now() - then))
        dp = times.toParts(since)
        prefix = strformat.fmt("{dp[Hours]}:{dp[Minutes]:02d}:{dp[Seconds]:02d}s ")
    write(stderr, prefix)
    log(words)

proc adjustThreadPool*(n: int) =
    ## n==0 => use ncpus
    ## n==-1 => do not alter threadpool size (to avoid a weird problem for now)
    log("(ThreadPool is currently not used.)")
    #var size = n
    #if n == 0:
    #    size = cpuinfo.countProcessors()
    #if size > threadpool.MaxThreadPoolSize:
    #    size = threadpool.MaxThreadPoolSize
    #if size == -1:
    #    log("ThreadPoolsize=", size,
    #        " (i.e. do not change)",
    #        ", MaxThreadPoolSize=", threadpool.MaxThreadPoolSize,
    #        ", NumCpus=", cpuinfo.countProcessors())
    #    return
    #log("ThreadPoolsize=", size,
    #    ", MaxThreadPoolSize=", threadpool.MaxThreadPoolSize,
    #    ", NumCpus=", cpuinfo.countProcessors())
    #threadpool.setMaxPoolSize(size)

iterator walk*(dir: string, followlinks = false, relative = false): string =
    ## similar to python os.walk(), but always topdown and no "onerror"
    # Slow! 30x slower than Unix find.
    let followFilter = if followLinks: {os.pcDir, os.pcLinkToDir} else: {os.pcDir}
    let yieldFilter = {os.pcFile, os.pcLinkToFile}
    for p in os.walkDirRec(dir, yieldFilter = yieldFilter,
            followFilter = followFilter, relative = relative):
        yield p

iterator readProc*(cmd: string): string =
    ## Stream from Unix subprocess, e.g. "find .".
    ## But if cmd=="-", stream directly from stdin.
    if cmd == "-":
        log("Reading from stdin...")
        for line in lines(stdin):
            yield line
    else:
        log("Reading from '" & cmd & "'...")
        var p = osproc.startProcess(cmd, options = {poEvalCommand})
        if osproc.peekExitCode(p) > 0:
            let msg = "Immedate failure in readProc startProcess('" & cmd & "')"
            raiseEx(msg)
        defer: osproc.close(p)
        for line in streams.lines(osproc.outputStream(p)):
            yield line

iterator readProcInMemory(cmd: string): string =
    ## Read from Unix subprocess, e.g. "find .", into memory.
    ## But if cmd=="-", stream directly from stdin.
    if cmd == "-":
        log("Reading from stdin...")
        for line in lines(stdin):
            yield line
    else:
        log("Reading from '" & cmd & "'...")
        let found = osproc.execProcess(cmd, options = {poEvalCommand})
        var sin = streams.newStringStream(found)
        for line in streams.lines(sin):
            yield line

proc removeFile*(fn: string, failIfMissing = false) =
    if failIfMissing and not os.existsFile(fn):
        raiseEx("Cannot remove non-existent file '" & fn & "'")
    log("rm -f ", fn)
    os.removeFile(fn)

proc removeFiles*(fns: openarray[string], failIfMissing = false) =
    for fn in fns:
        removeFile(fn, failIfMissing)

proc which*(exe: string) =
    let cmd = "which " & exe
    log(cmd)
    discard execCmd(cmd)

proc thousands*(v: SomeInteger): string =
    if v == 0:
        return "0"
    var i: type(v) = v
    let negative = (i < 0)
    i = abs(i)
    #result = strformat.fmt"{i mod 1000:03}"
    #i = i div 1000
    while i > 0:
        result = strformat.fmt"{i mod 1000:03}," & result
        i = i div 1000
    # Drop tailing comma.
    assert result[^1] == ','
    result = result[0 .. ^2]
    # Drop leading 0s.
    while result[0] == '0':
        result = result[1 .. ^1]
    if negative:
        result = '-' & result
