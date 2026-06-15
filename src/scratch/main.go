package main

import (
	"fmt"
)

func main() {

	bx := []byte("A01C46B31D11")

	fmt.Printf("Detected invalid barcode (%v), skipping.\n", string(bx))

}
