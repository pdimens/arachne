// Copyright © 2026 Pavel Dimens | Github: pdimens
package cmd

import (
	"arachne/aligner"
	"fmt"
	"os"
	"path/filepath"

	"github.com/spf13/cobra"
)

// alignCmd represents the align command
var alignCmd = &cobra.Command{
	Short:   "Align linked-read sequences to a reference",
	Use:     "align [flags] --sample-id REF.fa R1.fq R2.fq",
	Example: "arachne align -t 12 --sample-id sample1 ref.fa smp1.R1.fq.gz smp1.R2.fq.gz > smp1.sam",
	Long: "Align (short-read) linked-read sequences to a reference. Use \033[4;34marachne prep\033[0m to format " +
		"input FASTQ files for the aligner. Inputs must be sorted by barcode and in \"standard\" format, and can be Gzipped. " +
		"See the documentation for more information: https://pdimens.github.io/arachne",
	DisableFlagsInUseLine: true,
	SilenceUsage:          true,
	Args: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			fmt.Printf("%s", cmd.UsageString())
			return fmt.Errorf("please provide inputs")
		}
		if err := cobra.ExactArgs(3)(cmd, args); err != nil {
			return err
		}
		for _, j := range args {
			err := fileExists(j)
			if !err {
				return fmt.Errorf("file does not exist: %s", j)
			}
		}
		// if reference index files don't exist, run bwa index on reference
		exts := []string{".amb", ".ann", ".bwt", ".pac", ".sa"}
		for _, i := range exts {
			if _, err := os.Stat(args[0] + i); err != nil {
				return fmt.Errorf("Missing critical reference index file: %s (and possibly others). Please index reference with \033[94;1mbwa index\033[0m or (\033[94;1marachne index\033[0m", filepath.Base(args[0])+i)
			}
		}
		return nil
	},
	Run: arachneAlign,
}

func init() {
	rootCmd.AddCommand(alignCmd)

	//---Command line arguments-------------
	alignCmd.Flags().StringP("centromeres", "c", "", "TSV file describing known centromeres as CEN<chrname> <chrname> <start> <stop>")
	alignCmd.Flags().BoolP("comments", "C", false, "Append comments to SAM output")
	alignCmd.Flags().Float64P("improper-pair-penalty", "i", 4.0, "Penalty for improper pair")
	alignCmd.Flags().Int64P("infer-distance", "d", 50000, "Distance at which to consider reads with the same barcode to originate from different molecules")
	alignCmd.Flags().StringP("sample-id", "s", "", "Sample name (required)")
	if err := alignCmd.MarkFlagRequired("sample-id"); err != nil {
		panic(err)
	}
	alignCmd.Flags().IntP("threads", "t", 4, "Threads to use")
	alignCmd.Flags().BoolP("verbose", "v", false, "Verbose output")
}

func arachneAlign(cmd *cobra.Command, args []string) {
	//---Flag validations and failsafes -------------
	var debugSpoof bool
	sampleID, err := cmd.Flags().GetString("sample-id")
	if err != nil {
		panic(err)
	}

	inferDistance, err := cmd.Flags().GetInt64("infer-distance")
	if err != nil {
		panic(err)
	}
	inferDistance = max(inferDistance, 100)

	centromeres, err := cmd.Flags().GetString("centromeres")
	if err != nil {
		panic(err)
	}
	if centromeres != "" {
		err := fileExists(centromeres)
		if !err {
			panic(fmt.Errorf("file does not exist: %s", centromeres))
		}
	}

	threads, err := cmd.Flags().GetInt("threads")
	if err != nil {
		panic(err)
	}
	threads = max(threads, 1)

	improperPairPenalty, err := cmd.Flags().GetFloat64("improper-pair-penalty")
	if err != nil {
		panic(err)
	}
	if improperPairPenalty < 0.0 {
		improperPairPenalty *= -1.0
	}

	verbose, err := cmd.Flags().GetBool("verbose")
	if err != nil {
		panic(err)
	}
	//--- Setup config and run --------------------
	config := aligner.ArachneArgs{
		Reference:             &args[0],
		R1:                    &args[1],
		R2:                    &args[2],
		Improper_pair_penalty: &improperPairPenalty,
		Sample_id:             &sampleID,
		Threads:               &threads,
		DEBUG:                 &debugSpoof,
		InferDistance:         &inferDistance,
		DebugTags:             &debugSpoof,
		DebugPrintMove:        &debugSpoof,
		Centromeres:           &centromeres,
		Verbose:               &verbose,
	}
	aligner.Arachne(config)
}
