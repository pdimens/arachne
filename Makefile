export CGO_LDFLAGS = -L$(shell pwd)/gobwa/bwa

VERSION=0.2

all: arachne gobwa/bwa/libbwa.a gobwa/bwa/bwa

arachne: gobwa/bwa/libbwa.a gobwa/bwa/bwa
	@echo "Building arachne"
	mkdir -p bin/
	go build -ldflags "-X arachne/aligner.VERSION=$(VERSION) -s -w" -o bin/$@
	cp gobwa/bwa/bwa bin/
	chmod +x bin/arachne

gobwa/bwa/libbwa.a gobwa/bwa/bwa &:
	@echo "Building BWA"
	$(MAKE) -j 4 -C gobwa/bwa libbwa.a bwa

# --- experimental: cgo binding to lh3/minibwa (not part of `all`) ---------
gominibwa/minibwa/libminibwa.a gominibwa/minibwa/minibwa &:
	@echo "Building minibwa"
	$(MAKE) -j 4 -C gominibwa/minibwa libminibwa.a minibwa

test-gominibwa: export CGO_LDFLAGS = -L$(shell pwd)/gominibwa/minibwa
test-gominibwa: gominibwa/minibwa/libminibwa.a gominibwa/minibwa/minibwa
	go test -tags minibwa ./gominibwa

.PHONY: test-gominibwa

clean:
	@echo "Cleaning Build"
	rm -Rf bin/
	$(MAKE) -j 4 -C gobwa/bwa clean
	-$(MAKE) -C gominibwa/minibwa clean
