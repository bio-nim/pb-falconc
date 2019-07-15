# Falcon C Utilities
 TODO: WHAT do they do?

## Requirements
- **nim (0.20.0 or greater)**
- **htslib** at runtime
- **zlib** at runtime (via `dlopen()`)

## Set up
    source bash/module.sh
    export NIMBLE_DIR=$(pwd)/.git/NIMBLE_DIR  # e.g.

## Two universes

### Internet

With web access, we can use **nimble** to install **nim** dependencies.

However, one dependency is not yet in nimble: **networkx**.
You can "nimble install" networkx via

    make -f nimble.makefile sub

(We should probably copy networkx under src/ until it is a nimble-package. TODO.)

After that, all the standard "nimble" commands work. The `nimble.makefile`
shows how to run those, but `nimble --help` should suffice.

### Non-internet

When we build for mobs or bamboo, we lack web access.

All the **nim** dependencies are under `vendor/` as
"subtrees". (See `vendor/readme.md`)

To install from "vendor", use

    make rsync

Then, the standard "nimble" commands will *not* work,
but the standard "nim" commands will!

    make rsync
    nim c src/falconc.nim

Or simply:

    make test
    make integ  # currently broken
    make build
    make install

### Note: No mixing!

You must use *either* `nimble` *or* `rsync/nim`.

"nimble" will not work properly on the rsynced vendor packages.
