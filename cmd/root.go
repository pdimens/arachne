// Copyright © 2026 Pavel Dimens | Github: pdimens
package cmd

import (
	"arachne/aligner"
	"os"

	"github.com/spf13/cobra"
)

// rootCmd represents the base command when called without any subcommands
var rootCmd = &cobra.Command{
	Use:     "arachne",
	Version: aligner.VERSION,
	Short:   "A linked-read sequence aligner",
	Long: "Arachne linked-read aligner\n\nArachne is a successor to 10X Genomics' Lariat aligner that " +
		"is compatible with all current linked-read technologies (e.g. haplotagging, stLFR, TELLseq). " +
		"Input FASTQ files can come from any linked-read technology, provided they:\n" +
		"  - are a set of paired-end reads\n" +
		"  - are in \"standard\" format\n" +
		"    - have barcodes in a \033[92;1mBX:Z\033[0m SAM tag\n" +
		"      - e.g. \033[92;1mBX:Z:ATGGACTAGA\033[0m\n" +
		"    - have barcode validations (0|1) in a \033[92;1mVX:i\033[0m SAM tag\n" +
		"      - \033[92;1mVX:i:1\033[0m if valid, \033[92;1mVX:i:0\033[0m if invalid\n" +
		"  - are sorted by barcode\n\n",
}

func init() {
	rootCmd.CompletionOptions.DisableDefaultCmd = true
	rootCmd.SetHelpCommand(&cobra.Command{Hidden: true})
}

// Execute adds all child commands to the root command and sets flags appropriately.
// This is called by main.main(). It only needs to happen once to the rootCmd.
func Execute() {
	err := rootCmd.Execute()
	if err != nil {
		os.Exit(1)
	}
}
