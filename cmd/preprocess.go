// Copyright © 2026 Pavel Dimens | Github: pdimens
package cmd

import (
	"arachne/preprocess"
	"fmt"

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
	Args: func(cmd *cobra.Command, args []string) error {
		// Optionally run one of the validators provided by cobra
		if err := cobra.ExactArgs(3)(cmd, args); err != nil {
			return err
		}
		_, err := preprocess.FileExists(args[1])
		if err != nil {
			return err
		}
		_, err = preprocess.FileExists(args[2])
		if err != nil {
			return err
		}
		return nil
	},
	Run: func(cmd *cobra.Command, args []string) {
		threads, err := cmd.Flags().GetInt("threads")
		if err != nil {
			panic(err)
		}
		threads = max(threads, 1)

		fmt.Println("preprocess called with", threads, args[0], args[1], args[2])
	},
}

func init() {
	rootCmd.AddCommand(preprocessCmd)
	preprocessCmd.Flags().IntP("threads", "t", 4, "Threads to use")
}
