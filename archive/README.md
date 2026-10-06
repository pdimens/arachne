# archive

Retired code kept for reference. Nothing here is built, tested, or shipped.

## gobwa

The original cgo binding to [bwa](https://github.com/lh3/bwa) (`bwa mem` seeding, chaining and
mate rescue, with `bwa` as the `bwa/` submodule). Arachne used it as its aligner until it was
replaced by minibwa (see `../gominibwa`).

Its Go files carry the `bwa_archive` build tag, so `go build ./...` and `go test ./...` skip them.
Reviving it means more than adding the tag: `GetChains` and `GetAlignments` in `aligner/` were
written against this package's chain-then-Smith-Waterman API and were replaced by a minibwa
implementation, so they would have to be restored from git history (the commit that introduced
this directory) before the aligner could use it again.
