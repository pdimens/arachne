package preprocess

import "math/rand"

// Source - https://stackoverflow.com/a/31832326
// Posted by icza, modified by community. See post 'Timeline' for change history
// Retrieved 2026-06-25, License - CC BY-SA 4.0
const letterBytes = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"

// Create a random 4-letter alphanumeric string to use as a temporary filename for samtools sort
func RandString() string {
	b := make([]byte, 4)
	for i := range b {
		b[i] = letterBytes[rand.Intn(len(letterBytes))]
	}
	return string(b)
}
