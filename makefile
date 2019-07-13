PREFIX?=${CURDIR}
NIMBLE_DIR?=${CURDIR}/nimbleDir
export NIMBLE_DIR
# or use --nimbleDir:${NIMBLE_DIR} everywhere
NIMBLE_INSTALL=nimble install --debug -y

D=/pbi/flash/cdunn/bb/dipsim
run:
	./falconc phasr -a $D/alignments/aln.0.001_0.001.sort.bam -r $D/rangen/random.fa -o results

all:
	${MAKE} sub
	${NIMBLE_INSTALL}
quick:
	nim c -r tests/t_kmers.nim
integ:
	${MAKE} -C integ-tests
	nimble integ --debug # slow, for now
help:
	nimble -h
	nimble tasks
test:
	nimble test --debug # uses "tests/" directory by default
sub:
	#cd repos/cligen; nimble develop  # so we never need to reinstall after edits
	#cd vendor/nim-kmers; ${NIMBLE_INSTALL}
	cd vendor/nim-networkx; ${NIMBLE_INSTALL}
pretty:
	find src tests integ-tests -name '*.nim' | xargs -L1 nimpretty --indent=4
rsync:
	mkdir -p ${NIMBLE_DIR}/pkgs/
	rsync -av vendor/nim-networkx/src/ ${NIMBLE_DIR}/pkgs/networkx-1.0.0/
	rsync -av vendor/nim-heap/ ${NIMBLE_DIR}/pkgs/binaryheap-0.1.1/
	rsync -av vendor/hts-nim/src/ ${NIMBLE_DIR}/pkgs/hts-0.2.15/
	rsync -av vendor/cligen/ ${NIMBLE_DIR}/pkgs/cligen-0.9.34/

# These 3 rules are for mobs/bamboo:
# Someday maybe --nimcache:${CURDIR}/.cache-nim
build:
	# We need a no-internet flag for "nimble install".
	# For now, we rsync and install manually.
	${MAKE} rsync
	nim c --listCmd -d:release src/falconc.nim # uses NIMBLE_DIR
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
