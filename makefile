PREFIX?=${CURDIR}
NIMBLE_DIR?=${CURDIR}/nimbleDir
export NIMBLE_DIR
# or use --nimbleDir:${NIMBLE_DIR} everywhere
NIMBLE_INSTALL=nimble install --debug -y

default:
	${NIMBLE_INSTALL}
test:
	nimble test --debug # uses "tests/" directory by default
pretty:
	find . -name '*.nim' | xargs -L1 nimpretty --indent=4

# These 3 rules are for mobs:
# Someday maybe --nimcache:${CURDIR}/.cache-nim
build:
	mkdir -p ${NIMBLE_DIR}/pkgs/
	# We need a no-internet flag for "nimble install".
	# For now, we install manually.
	rsync -av vendor/nim-networkx/src/ ${NIMBLE_DIR}/pkgs/networkx-1.0.0/
	rsync -av vendor/nim-heap/ ${NIMBLE_DIR}/pkgs/binaryheap-0.1.1/
	rsync -av vendor/hts-nim/src/ ${NIMBLE_DIR}/pkgs/hts-0.2.15/
	rsync -av vendor/cligen/ ${NIMBLE_DIR}/pkgs/cligen-0.9.34/
	nim c --listCmd src/falconc.nim # uses NIMBLE_DIR
install:
	mkdir -p ${PREFIX}/bin
	mv -f src/falconc ${PREFIX}/bin
clean:
	rm -rf ~/.cache/nim/falconc* # race-condition here
	rm -rf src/falconc ${NIMBLE_DIR}
