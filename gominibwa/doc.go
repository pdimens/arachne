// Package gominibwa is a cgo binding to lh3/minibwa
// (https://github.com/lh3/minibwa), the counterpart of package gobwa for
// bwa. It is experimental and not yet used by the aligner.
//
// Unlike bwa, where arachne first obtains chains from the library and then
// runs Smith-Waterman on them separately, minibwa returns complete
// alignments (with CIGARs) from a single call, so the API here is
// correspondingly smaller: load an Index, create a Mapper per goroutine,
// and call Map, MapPair or MapPairs.
//
// The binding is only compiled with the "minibwa" build tag so that the rest
// of arachne builds without libminibwa:
//
//	make -C gominibwa/minibwa libminibwa.a minibwa
//	CGO_LDFLAGS="-L$PWD/gominibwa/minibwa" go test -tags minibwa ./gominibwa
//
// Indexes are built with the minibwa binary (see Build) and consist of
// <ref>.l2b and <ref>.mbw.
package gominibwa
