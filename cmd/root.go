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
		"is compatible with all current linked-read technologies (e.g. haplotagging, stLFR, TELLseq), provided they are preprocessed correctly.",
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
