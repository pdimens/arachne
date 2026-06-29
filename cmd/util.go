package cmd

import (
	"fmt"
	"os"
	"os/exec"
)

func filecheck(filename string) error {
	//var err error
	if err := fileExists(filename); err != nil {
		return err
	}
	if err := fileReadable(filename); err != nil {
		return err
	}
	return nil
}

// check if file exists. Returns a true if it does, otherwise false.
func fileExists(filename string) error {
	info, err := os.Stat(filename)
	if os.IsNotExist(err) {
		return fmt.Errorf("\033[33;1m%s\033[0m does not exist", filename)
	}
	if info.IsDir() {
		return fmt.Errorf("\033[33;1m%s\033[0m is a directory", filename)
	}
	return nil
}

func fileReadable(filename string) error {
	file, err := os.Open(filename)
	if err != nil {
		return fmt.Errorf("\033[33;1m%s\033[0m does not have read persmissions.\n", filename)
	}
	file.Close()
	return nil
}

// check if executable 'e' is in path, returns true if it is, otherwise false
func checkIfExecInPath(e string) bool {
	// ignore output, only use error
	_, err := exec.LookPath(e)
	// return true if error nil, else false
	return err == nil
}
