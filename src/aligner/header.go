package aligner

import (
	"arachne/src/gobwa"
	"log"
	"os"
	"strings"
	"time"

	"github.com/biogo/hts/sam"
)

// Return a SAM header built from the reference, RG, sample ID, and arachne version
func buildHeader(ref *gobwa.GoBwaReference, read_groups, sample_id, version string) (*sam.Header, map[string]*sam.Reference) {
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

	for _, rg_id := range strings.Split(read_groups, ",") {
		// currently, the ID is composed of:
		// sample:library:gem_group:flowcell:lane
		rg_fields := strings.Split(rg_id, ":")
		if len(rg_fields) == 0 {
			log.Printf("Empty RG was specified, skipping")
		} else if len(rg_fields) < 5 {
			log.Printf("RG is not fully specified, skipping: %s", rg_id)
		} else {
			rg, err := sam.NewReadGroup(
				rg_id,                         //ID
				"",                            //CN
				"",                            //DS
				rg_fields[1]+"."+rg_fields[2], //LB = (input library).(gem group)
				"",                            //PG
				"ILLUMINA",                    //PL
				rg_id,                         //PU: just make same as ID?
				rg_fields[0],                  //SM
				"",
				"",
				time.Now(),
				0)
			if err != nil {
				panic(err)
			}
			h.AddReadGroup(rg)
		}
	}

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
