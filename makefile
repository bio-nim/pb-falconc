PREFIX?=${CURDIR}
NIMBLE_DIR?=${CURDIR}/nimbleDir
export NIMBLE_DIR
# or use --nimbleDir:${NIMBLE_DIR} everywhere
NIMBLE_INSTALL=nimble install --debug -y

install:
	${NIMBLE_INSTALL}
test:
	nimble test --debug # uses "tests/" directory by default
pretty:
	find . -name '*.nim' | xargs -L1 nimpretty --indent=4
