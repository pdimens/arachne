package preprocess

import (
	"arachne/fastqreader"
	"flag"
	"fmt"
	"io"
	"os"

	"github.com/shenwei356/xopen"
)

func main() {
	var threads int
	var activeRecord fastqreader.FastQRecord

	//----Command line arguments-----------------------
	flag.IntVar(&threads, "threads", 8, "Number of threads")
	flag.IntVar(&threads, "t", 8, "Number of threads")

	flag.Usage = func() {
		fmt.Fprint(os.Stderr, "Preprocess FASTQ file pair for Arachne linked-read aligner\n")
		fmt.Fprint(os.Stderr, "\n\033[94;1mUsage:\033[0m arachne-pre PREFIX sample.R1.fq sample.R2.fq\n")
	}

	prefix := flag.Arg(0)

	r1 := flag.Arg(1)
	fastqreader.FileExists(r1, "FASTQ")

	r2 := flag.Arg(2)
	fastqreader.FileExists(r2, "FASTQ")

	fastq, err := fastqreader.OpenFastQPair(r1, r2)
	if err != nil {
		panic(err)
	}
	defer fastq.Close()
	outfq1, err := xopen.Wopen(prefix + ".R1.fq.gz")
	if err != nil {
		panic(err)
	}
	outfq2, err := xopen.Wopen(prefix + ".R2.fq.gz")
	if err != nil {
		panic(err)
	}
	filtfq1, err := xopen.Wopen(prefix + ".removed.R1.fq.gz")
	if err != nil {
		panic(err)
	}
	filtfq2, err := xopen.Wopen(prefix + ".removed.R2.fq.gz")
	if err != nil {
		panic(err)
	}

	// ── loop through records ─────────────────────────────────────────────────────────
	for {
		err = fastq.ReadOneRecord(&activeRecord)
		if err != nil {
			// no more records
			if err == io.EOF {
				break
			} else {
				panic(err)
			}
		}

	}

}
