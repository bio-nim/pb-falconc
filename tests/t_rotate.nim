# vim: sw=4 ts=4 sts=4 tw=0 et:

import falconcpkg/rotate
import os
#import sets
import unittest

template withFile(f: untyped, filename: string, mode: FileMode,
                  body: untyped) =
    # from https://nim-lang.org/docs/tut2.html
    let fn = filename
    var f: File
    if open(f, fn, mode):
        try:
            body
        finally:
            close(f)
    else:
        quit("cannot open '" & fn & "' in mode " & $mode)


suite "rotate":
    #unittest.abortOnError = true # so we can examine files on error
    let subdir = "rotate-testdir"
    echo "CREATE ", subdir
    discard os.existsOrCreateDir(subdir)
    os.setCurrentDir(subdir)
    test "loadWhiteList missing":
        let wl = loadWhiteList("")
        check whitelisted(wl, "foo")
        check whitelisted(wl, "")
    test "loadWhiteList empty":
        let fn = "empty.txt"
        withFile(fout, fn, fmWrite):
            fout.write() # just touch
        let wl = loadWhiteList(fn)
        check not whitelisted(wl, "foo")
        check not whitelisted(wl, "")
    test "loadWhiteList 1-line":
        let fn = "sumpinsumpin.txt"
        withFile(fout, fn, fmWrite):
            fout.write("foo")
        let wl = loadWhiteList(fn)
        check whitelisted(wl, "foo")
        check not whitelisted(wl, "bar")
        check not whitelisted(wl, "")

    # Might skip all this if abortOnError.
    os.setCurrentDir("..")
    echo "REMOVE ", subdir
    os.removeDir(subdir)

suite "rotate-misc":
    test "chromFromHeader":
        assert chromFromHeader("@foo") == "foo"
        assert chromFromHeader(">foo") == "foo"
        assert chromFromHeader("@foo bar") == "foo"
        assert chromFromHeader(">foo bar") == "foo"
