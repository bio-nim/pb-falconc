# Changelog - IPA assembler

## Active version in development
### Changes
- No changes

## v1.13.0 - SL - Release 10.1.0
### Version
- `falconc` - `commit 235b6add85aea15c9de6e30c207e7a8329cd1c88 (origin/release/prep)` (Dec 16, 2020), `falconc version=1.13.0+git., nim-version=1.2.0`

### Changes
- Added a PAF splitter and indexer to be used with the new mapping via read tracking
- PAF splitter for the read-to-contig tracking-based-mapping for polishing.
- Subtool to rename IPA contig headers to the Falcon-compatible format, so that the contigs can be used with Falcon-Phase.
- Exposed the min average coverage cutoff which was previously hardcoded to 5. (For `m4filt*` subtools; new parameter `--min-avg-cov`.)
- New subtool `gfftools-gffsubtract`.
- Updated the `hts-nim` vendor library.

## v1.10.1 - SL - Release 10.0.0 (SL-release-10.0.0)
### Version
- `falconc` - commit 16390377ed98c528792c0ee405d81538913ad074 (tag: SL-release-10.0.0, origin/release/10.0.0) (Sep 23, 2020), `falconc version=1.10.1+git., nim-version=1.2.0`
