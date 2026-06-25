// Copyright © 2026 Pavel Dimens | Github: pdimens
package cmd

import (
	"fmt"
	"runtime"

	"arachne/preprocess"

	"github.com/spf13/cobra"
)

// preprocessCmd represents the preprocess command
var preprocessCmd = &cobra.Command{
	Use:     "preprocess [-t] PREFIX R1.fq R2.fq",
	Short:   "Format FASTQ files for alignment \033[94;1m[start here]\033[0m",
	Example: "preprocess -t 12 sample1 smp1.R1.fq.gz smp1.R2.fq.gz",
	Long: "Converts a set of paired-end FASTQ files into the format required for the arachne aligner. " +
		"For arachne to work correctly, input FASTQ files need to be properly paired, in \"standard\" format " +
		"(\033[94;1mBX:Z\033[0m and \033[94;1mVX:i\033[0m tags), and sorted by barcode. " +
		"Reads with invalid barcodes will be preserved separately so they can be aligned using another tool like BWA. " +
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
		maxCores := runtime.NumCPU()
		// clamp between 1 and max system threads
		threads = min(maxCores, max(threads, 1))
		runtime.GOMAXPROCS(threads)
		preprocess.Preprocess(threads, args[0], args[1], args[2])
		return nil
	},
}

func init() {
	rootCmd.AddCommand(preprocessCmd)
	preprocessCmd.Flags().IntP("threads", "t", 2, "Number of threads to use")
}
