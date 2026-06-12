package aligner

import (
	"arachne/src/gobwa"
	"fmt"
	"log"
	"os"
)

func bwaindex(fasta string) {
	fmt.Printf("Building BWA index for %s...\n", fasta)
	err := gobwa.GoBwaIndex(fasta)
	if err != nil {
		log.Fatalf("Error: Failed to build BWA reference index.\n%v", err)
		os.Exit(1)
	} else {
		fmt.Println("Index built successfully!")
	}
}
