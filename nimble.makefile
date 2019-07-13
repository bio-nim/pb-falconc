include common.makefile

PREFIX?=${CURDIR}
#NIMBLE_DIR?=${CURDIR}/nimbleDir
#export NIMBLE_DIR
# or use --nimbleDir:${NIMBLE_DIR} everywhere
NIMBLE_INSTALL=nimble install --debug -y

all:
	${MAKE} sub
	${NIMBLE_INSTALL}
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

.PHONY: test
