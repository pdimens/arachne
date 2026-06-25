package aligner

import (
	"os"
	"strings"
	"time"

	"arachne/gobwa"

	"github.com/biogo/hts/sam"
)

// Return a SAM header built from the reference, sample ID, and arachne version
func buildHeader(ref *gobwa.GoBwaReference, sampleid, version string) (*sam.Header, map[string]*sam.Reference) {
	contigs := make(map[string]*sam.Reference)
	references := make([]*sam.Reference, 0)

	gobwa.EnumerateContigs(ref, func(name string, length int) {
		r, err := sam.NewReference(name, name, "NA", length, nil, nil)
		if err != nil {
			panic(err)
		}
		references = append(references, r)
		contigs[name] = r
	})

	h, err := sam.NewHeader([]byte(""), references)
	if err != nil {
		panic(err)
	}

	rg, err := sam.NewReadGroup(
		sampleid,   //ID
		"",         //CN
		"",         //DS
		"",         //LB
		"",         //PG
		"ILLUMINA", //PL
		"",         //PU: just make same as ID?
		sampleid,   //SM
		"",
		"",
		time.Now(),
		0)
	if err != nil {
		panic(err)
	}
	h.AddReadGroup(rg)

	// Add a program line for arachne
	prog := sam.NewProgram(
		"arachne",                  // ID
		"arachne",                  // PN
		strings.Join(os.Args, " "), // CL
		"",                         // PP - empty bc Arachne produces the initial SAM
		version)                    // VN
	h.AddProgram(prog)
	return h, contigs
}
