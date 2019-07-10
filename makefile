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
build:
	cd vendor/nim-networkx; nimble install
	nimble build
install:
	mkdir -p ${PREFIX}/bin
	mv -f ./falconc ${PREFIX}/bin
clean:
	rm -rf ./falconc ${NIMBLE_DIR}
