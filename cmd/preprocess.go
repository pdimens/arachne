// Copyright © 2026 Pavel Dimens | Github: pdimens
package cmd

import (
	"fmt"
	"runtime"

	"arachne/preprocess"

	"github.com/spf13/cobra"
)

// preCmd represents the preprocess command
var prep = &cobra.Command{
	Use:     "prep [-t] PREFIX R1.fq R2.fq",
	Short:   "Format FASTQ files for alignment \033[94;1m[start here]\033[0m",
	Example: "prep -t 12 sample1 smp1.R1.fq.gz smp1.R2.fq.gz",
	Long: "Sorts a pair of \"standard\"-format (\033[94;1mBX:Z\033[0m and \033[94;1mVX:i\033[0m tags) linked-read FASTQ files " +
		"by barcode, as required by the arachne aligner. Input must already be in standard format; " +
		"use \033[94;1mdjinn\033[0m to convert haplotagging, stLFR, or TELLseq FASTQ data first. " +
		"\033[4;32mRequires samtools to be available in your PATH.\033[0m",
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
		if err := filecheck(args[1]); err != nil {
			return err
		}
		if err := filecheck(args[2]); err != nil {
			return err
		}
		if !checkIfExecInPath("samtools") {
			return fmt.Errorf("samtools was not found on the PATH:")
		}
		return nil
	},
	RunE: func(cmd *cobra.Command, args []string) error {
		threads, err := cmd.Flags().GetInt("threads")
		if err != nil {
			return err
		}
		maxCores := runtime.NumCPU()
		// clamp between 1 and max system threads
		threads = min(maxCores, max(threads, 1))
		runtime.GOMAXPROCS(threads)
		return preprocess.Preprocess(threads, args[0], args[1], args[2])
	},
}

func init() {
	rootCmd.AddCommand(prep)
	prep.Flags().IntP("threads", "t", 2, "Number of threads to use")
}
