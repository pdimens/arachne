package cmd

import (
	"arachne/gobwa"
	"fmt"

	"github.com/spf13/cobra"
)

var indexCmd = &cobra.Command{
	Short:                 "Index reference fasta",
	Use:                   "index REF.fa",
	Example:               "arachne index ref.fa",
	Long:                  "This is a small wrapper to use BWA (included) to index a reference FASTA file prior to alignment.",
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
		if !fileExists(args[0]) {
			return fmt.Errorf("file does not exist: %s", args[0])
		}
		return nil
	},
	RunE: func(cmd *cobra.Command, args []string) error {
		return bwaindex(args[0])
	},
}

func init() {
	rootCmd.AddCommand(indexCmd)
}

// Thin shim/wrapper to index a reference with bwa index. Stderr only prints
// on error.
func bwaindex(fasta string) error {
	fmt.Printf("   Building BWA index for %s...\n", fasta)
	err := gobwa.GoBwaIndex(fasta)
	if err != nil {
		return fmt.Errorf("Error: Failed to build BWA reference index.\n%v", err)
	}
	fmt.Println("   Index built successfully!")
	return nil
}
