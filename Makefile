export CGO_LDFLAGS = -L$(shell pwd)/gobwa/bwa -L$(shell pwd)/jemalloc
export GOPATH=$(shell pwd)

VERSION=1.0-dev

GO_VERSION=$(strip $(shell go version | sed 's/.*go\([0-9]*\.[0-9]*\).*/\1/'))

all: arachne gobwa/bwa/libbwa.a gobwa/bwa/bwa

gobwa/bwa/bwa: gobwa/bwa/libbwa.a
	@echo "Building bwa binary (for bwa index)"
	make -C gobwa/bwa bwa

arachne: gobwa/bwa/libbwa.a gobwa/bwa/bwa
	@echo "Building arachne"
	mkdir -p bin/
	go build -o bin/arachne $@
	cp gobwa/bwa/bwa bin/
	chmod +x bin/arachne

gobwa/bwa/libbwa.a:
	@echo "Building BWA"
	make -C gobwa/bwa libbwa.a

clean:
	@echo "Cleaning Build"
	rm -Rf bin/
	$(MAKE) -C gobwa/bwa clean
