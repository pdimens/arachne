export CGO_LDFLAGS = -L$(shell pwd)/src/gobwa/bwa -L$(shell pwd)/src/jemalloc
export GOPATH=$(shell pwd)

VERSION=1.0-dev

GO_VERSION=$(strip $(shell go version | sed 's/.*go\([0-9]*\.[0-9]*\).*/\1/'))

all: arachne src/gobwa/bwa/libbwa.a src/gobwa/bwa/bwa

src/gobwa/bwa/bwa: src/gobwa/bwa/libbwa.a
	@echo "Building bwa binary (for bwa index)"
	make -C src/gobwa/bwa bwa

arachne: src/gobwa/bwa/libbwa.a src/gobwa/bwa/bwa
	@echo "Building arachne"
	mkdir -p bin/
	go build -o bin/arachne $@
	cp src/gobwa/bwa/bwa bin/
	chmod +x bin/arachne

src/gobwa/bwa/libbwa.a:
	@echo "Building BWA"
	make -C src/gobwa/bwa libbwa.a

clean:
	@echo "Cleaning Build"
	rm -Rf bin/
	$(MAKE) -C src/gobwa/bwa clean

test:
	cd src/test; go test -v
