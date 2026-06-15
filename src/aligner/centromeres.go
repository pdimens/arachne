package aligner

import (
	"bufio"
	"log"
	"os"
	"strconv"
	"strings"
)

type Region struct {
	start int
	end   int
}

func loadCentromeres(filename *string) map[string]Region {
	toRet := map[string]Region{}
	if *filename == "" {
		return toRet
	}
	file, _ := os.Open(*filename)
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, "CEN") {
			tokens := strings.Split(line, "\t")
			if len(tokens) < 4 {
				continue
			}
			chrom := tokens[1]
			start, err := strconv.Atoi(tokens[2])
			if err != nil {
				continue
			}
			end, err := strconv.Atoi(tokens[3])
			if err != nil {
				continue
			}
			toRet[chrom] = Region{start: start, end: end}
		}
	}
	if err := scanner.Err(); err != nil {
		log.Fatalf("Unable to read Centromeres file.\n%v", err)
	}
	return toRet
}
