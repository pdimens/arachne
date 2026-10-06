// Copyright © 2026 Pavel Dimens | Github: pdimens
package cmd

import (
	"arachne/gominibwa"
	"fmt"

	"github.com/spf13/cobra"
)

var indexCmd = &cobra.Command{
	Short:   "Index reference.fasta",
	Use:     "index REF.fa",
	Example: "arachne index ref.fa",
	Long: "A small wrapper to use minibwa (included) to index a reference FASTA file prior to alignment. " +
		"This creates REF.fa.l2b and REF.fa.mbw alongside the FASTA.",
	DisableFlagsInUseLine: true,
	SilenceUsage:          true,
	Args: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			fmt.Printf("%s", cmd.UsageString())
			return fmt.Errorf("please provide input reference FASTA")
		}
		if err := cobra.ExactArgs(1)(cmd, args); err != nil {
			return err
		}
		if err := filecheck(args[0]); err != nil {
			return err
		}
		return nil
	},
	RunE: func(cmd *cobra.Command, args []string) error {
		threads, _ := cmd.Flags().GetInt("threads")
		return minibwaIndex(args[0], max(threads, 1))
	},
}

func init() {
	rootCmd.AddCommand(indexCmd)
	indexCmd.Flags().IntP("threads", "t", 4, "Threads to use")
}

// Thin shim/wrapper to index a reference with minibwa index. Stderr only
// prints on error.
func minibwaIndex(fasta string, threads int) error {
	fmt.Printf("   Building minibwa index for %s...\n", fasta)
	err := gominibwa.Build(fasta, threads)
	if err != nil {
		return fmt.Errorf("Error: Failed to build minibwa reference index.\n%v", err)
	}
	fmt.Println("   Index built successfully!")
	return nil
}
