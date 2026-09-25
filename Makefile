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

clean:
	@echo "Cleaning Build"
	rm -Rf bin/
	$(MAKE) -j 4 -C gobwa/bwa clean
