// Copyright © 2026 Pavel Dimens | Github: pdimens
package cmd

import (
	"arachne/aligner"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"time"

	"github.com/spf13/cobra"
)

// alignCmd represents the align command
var alignCmd = &cobra.Command{
	Short:   "Align linked-read sequences to a reference",
	Use:     "align [flags] --sample-id REF.fa R1.fq R2.fq",
	Example: "arachne align -t 12 --sample-id sample1 ref.fa smp1.R1.fq.gz smp1.R2.fq.gz > smp1.bam",
	Long: "Align (short-read) linked-read sequences to a reference. " +
		"Inputs must be in 'standard' format (use \033[4;34mdjinn\033[0m) and sorted by barcode (use \033[4;34marachne prep\033[0m). " +
		"Output is uncompressed BAM.\n\n" +
		"A --centromeres file is tab-delimited BED format <chrname> <start> <stop>. Sequences that align to centromeric regions will have " +
		"their MAPQ set to 0, as mapping to centromeric regions is unreliable.\n" +
		"Documentation: https://pdimens.github.io/arachne",
	DisableFlagsInUseLine: true,
	SilenceUsage:          true,
	Args: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			fmt.Printf("%s", cmd.UsageString())
		}
		if err := cobra.ExactArgs(3)(cmd, args); err != nil {
			return err
		}
		for _, j := range args {
			if err := filecheck(j); err != nil {
				return err
			}
		}
		// if reference index files don't exist, run bwa index on reference
		exts := []string{".amb", ".ann", ".bwt", ".pac", ".sa"}
		for _, i := range exts {
			if _, err := os.Stat(args[0] + i); err != nil {
				return fmt.Errorf("missing reference index file: %s\nPlease index reference with \033[94;1marachne index\033[0m or (\033[94;1mbwa index\033[0m)", filepath.Base(args[0])+i)
			}
		}
		return nil
	},
	RunE: arachneAlign,
}

// normalizeImproperPairPenalty ensures the improper-pair penalty is always
// applied as a penalty, never a bonus. aligner.scoreAlignment() *adds* this
// value to the pair score, so it must be <= 0 regardless of the sign the
// user passed on the command line.
func normalizeImproperPairPenalty(v float64) float64 {
	return -math.Abs(v)
}

func init() {
	rootCmd.AddCommand(alignCmd)

	//---Command line arguments-------------
	alignCmd.Flags().StringP("centromeres", "c", "", "BED file describing known centromeres (optional, see --help)")
	alignCmd.Flags().BoolP("comments", "C", false, "Append comments (non-BX/VX) to SAM output")
	alignCmd.Flags().Float64P("improper-pair-penalty", "i", 4.0, "Penalty for improper pair")
	alignCmd.Flags().Int64P("infer-distance", "d", 50000, "Distance at which to consider reads with the same barcode to be from different molecules")
	alignCmd.Flags().StringP("sample-id", "s", "", "Sample name (required)")
	if err := alignCmd.MarkFlagRequired("sample-id"); err != nil {
		panic(err)
	}
	alignCmd.Flags().IntP("threads", "t", 4, "Threads to use")
	alignCmd.Flags().BoolP("verbose", "v", false, "Verbose output")
}

func arachneAlign(cmd *cobra.Command, args []string) error {
	//---Flag validations and failsafes -------------
	var debugSpoof bool
	sampleID, err := cmd.Flags().GetString("sample-id")
	if err != nil {
		return err
	}

	inferDistance, err := cmd.Flags().GetInt64("infer-distance")
	if err != nil {
		return err
	}
	inferDistance = max(inferDistance, 100)

	centromeres, err := cmd.Flags().GetString("centromeres")
	if err != nil {
		return err
	}
	if centromeres != "" {
		if err := filecheck(centromeres); err != nil {
			return err
		}
	}

	threads, err := cmd.Flags().GetInt("threads")
	if err != nil {
		return err
	}
	threads = max(threads, 1)

	improperPairPenalty, err := cmd.Flags().GetFloat64("improper-pair-penalty")
	if err != nil {
		return err
	}
	improperPairPenalty = normalizeImproperPairPenalty(improperPairPenalty)

	comments, err := cmd.Flags().GetBool("comments")
	if err != nil {
		return err
	}

	verbose, err := cmd.Flags().GetBool("verbose")
	if err != nil {
		return err
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
		Comments:              &comments,
	}
	start := time.Now()
	fmt.Fprintf(os.Stderr, "🕷️  Starting arachne. Version: %s\n", aligner.VERSION)

	aligner.Arachne(config)

	elapsed := time.Since(start).Round(time.Second).String()
	fmt.Fprintf(os.Stderr, "🕸️  Arachne finished successfully! Elapsed: %s\n\n", elapsed)
	return nil
}
