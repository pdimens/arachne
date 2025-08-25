export CGO_LDFLAGS = -L$(shell pwd)/src/gobwa/bwa -L$(shell pwd)/src/jemalloc
export GOPATH=$(shell pwd)

VERSION=$(shell git describe --tags --always --dirty)

GO_VERSION=$(strip $(shell go version | sed 's/.*go\([0-9]*\.[0-9]*\).*/\1/'))

all: arachne arachne-preprocess

arachne-preprocess:
	@echo "Building arachne-preprocess"
	go install ./src/preprocess

arachne: src/gobwa/bwa/libbwa.a
	@echo "Building arachne"
	go install -ldflags "-X aligner.__VERSION__='$(VERSION)'" $@

src/gobwa/bwa/libbwa.a:
	make -C src/gobwa/bwa libbwa.a

clean:
	rm -Rf bin/ pkg
	$(MAKE) -C src/gobwa/bwa clean

test:
	cd src/test; go test -v