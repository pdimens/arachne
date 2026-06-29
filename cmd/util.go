package cmd

import (
	"os"
	"os/exec"
)

// check if file exists. Returns a true if it does, otherwise false.
func fileExists(filename string) bool {
	info, err := os.Stat(filename)
	if os.IsNotExist(err) {
		return false
	}
	return !info.IsDir()
}

// check if executable 'e' is in path, returns true if it is, otherwise false
func checkIfExecInPath(e string) bool {
	// ignore output, only use error
	_, err := exec.LookPath(e)
	// return true if error nil, else false
	return err == nil
}
