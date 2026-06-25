// Copyright © 2026 Pavel Dimens | Github: pdimens
package cmd

import (
	"fmt"

	"arachne/preprocess"

	"github.com/spf13/cobra"
)

// preprocessCmd represents the preprocess command
var preprocessCmd = &cobra.Command{
	Use:     "preprocess [-t] PREFIX R1.fq R2.fq",
	Short:   "Format FASTQ files for alignment",
	Example: "preprocess -t 12 sample1 smp1.R1.fq.gz smp1.R2.fq.gz",
	Long: "The command converts a set of paired-end FASTQ files into the format required for the arachne aligner. " +
		"For arachne to work correctly, input FASTQ files need to be in \"standard\" format and sorted by barcode (at minimum). " +
		"Reads with invalid and singleton barcodes will be preserved separately so they can be aligned using another tool like BWA.",
	DisableFlagsInUseLine: true,
	SilenceUsage:          true,
	Args: func(cmd *cobra.Command, args []string) error {
		if err := cobra.ExactArgs(3)(cmd, args); err != nil {
			return err
		}
		err := fileExists(args[1])
		if !err {
			return fmt.Errorf("file does not exist: %s", args[1])
		}
		err = fileExists(args[2])
		if !err {
			return fmt.Errorf("file does not exist: %s", args[2])
		}
		return nil
	},
	RunE: func(cmd *cobra.Command, args []string) error {
		threads, err := cmd.Flags().GetInt("threads")
		if err != nil {
			return err
		}
		threads = max(threads, 1)
		preprocess.Preprocess(threads, args[0], args[1], args[2])
		return nil
	},
}

func init() {
	rootCmd.AddCommand(preprocessCmd)
	preprocessCmd.Flags().IntP("threads", "t", 2, fmt.Sprintf("Number of threads to use (min: %d)", 2))
}
