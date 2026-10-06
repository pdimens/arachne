// Package gominibwa is a cgo binding to lh3/minibwa
// (https://github.com/lh3/minibwa), arachne's aligner. It replaced the bwa
// binding, which is kept in archive/gobwa.
//
// Unlike bwa, where arachne used to obtain chains from the library and then
// run Smith-Waterman on them separately, minibwa returns complete
// alignments (with CIGARs) from a single call. The API is correspondingly
// small: load an Index, create a Mapper per goroutine, and call Map,
// MapPair or MapPairs.
//
// The minibwa submodule must be built before this package links:
//
//	make -C gominibwa/minibwa libminibwa.a minibwa
//	CGO_LDFLAGS="-L$PWD/gominibwa/minibwa" go test ./gominibwa
//
// (`make` at the repository root does this.) Indexes are built with the
// minibwa binary (see Build) and consist of <ref>.l2b and <ref>.mbw.
package gominibwa
