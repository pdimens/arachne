export CGO_LDFLAGS = -L$(shell pwd)/src/gobwa/bwa -L$(shell pwd)/src/jemalloc
export GOPATH=$(shell pwd)

VERSION=1.0-dev

GO_VERSION=$(strip $(shell go version | sed 's/.*go\([0-9]*\.[0-9]*\).*/\1/'))

all: arachne src/gobwa/bwa/libbwa.a

arachne: src/gobwa/bwa/libbwa.a
	@echo "Building arachne"
	mkdir -p bin/
	go build -ldflags "-X aligner.__VERSION__='$(VERSION)'" -o bin/arachne $@

src/gobwa/bwa/libbwa.a:
	@echo "Building BWA"
	make -C src/gobwa/bwa libbwa.a

clean:
	@echo "Cleaning Build"
	rm -Rf bin/ pkg
	$(MAKE) -C src/gobwa/bwa clean

test:
	cd src/test; go test -v
