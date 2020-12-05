include common.makefile

PREFIX?=${CURDIR}
NIMBLE_DIR?=${CURDIR}/nimbleDir
export NIMBLE_DIR
# or use --nimbleDir:${NIMBLE_DIR} everywhere
NIMBLE_INSTALL=nimble install --debug -y

default: build
nim:
	nim c --listCmd -d:release src/falconc.nim # uses NIMBLE_DIR
test:
	which nim
	nim -v
	${MAKE} -C tests/
integ:
	${MAKE} -C integ-tests/
rsync:
	mkdir -p ${NIMBLE_DIR}/pkgs/
	rsync -av vendor/nim-networkx/src/ ${NIMBLE_DIR}/pkgs/networkx-1.0.0/
	rsync -av vendor/nim-heap/ ${NIMBLE_DIR}/pkgs/binaryheap-0.1.1/
	rsync -av vendor/hts-nim/src/ ${NIMBLE_DIR}/pkgs/hts-0.3.14/
	rsync -av vendor/msgpack4nim/ ${NIMBLE_DIR}/pkgs/msgpack4nim-0.2.9/
	rsync -av vendor/cligen/ ${NIMBLE_DIR}/pkgs/cligen-1.3.2/
	rsync -av vendor/threadpools/ ${NIMBLE_DIR}/pkgs/threadpools-0.1.0/
	rsync -av vendor/comprehension/ ${NIMBLE_DIR}/pkgs/comprehension-0.1.0/
	rsync -av vendor/c_alikes/ ${NIMBLE_DIR}/pkgs/c_alikes-0.2.0/
	rsync -av vendor/BitVector/ ${NIMBLE_DIR}/pkgs/BitVector-0.4.10/

# These 3 rules are for mobs/bamboo:
# Someday maybe --nimcache:${CURDIR}/.cache-nim
#build: R=$(shell git rev-parse HEAD)
build:
	# We need a no-internet flag for "nimble install".
	# For now, we rsync and install manually.
	${MAKE} rsync
	# "nim c" uses NIMBLE_DIR
	nim c --listCmd --threads:on -d:release --styleCheck:hint -d:cGitSha1=$R src/falconc.nim
	#nim c --listCmd --threads:on -d:release --styleCheck:hint src/falconc.nim
install:
	mkdir -p ${PREFIX}/bin
	mv -f src/falconc ${PREFIX}/bin
clean:
	#rm -rf ~/.cache/nim/falconc* # race-condition here
	rm -rf src/falconc ${NIMBLE_DIR}

#These might control the cache-dir:
# XDG_CACHE_HOME
# XDG_CONFIG_HOME

.PHONY: test
