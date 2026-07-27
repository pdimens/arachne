export CGO_LDFLAGS = -L$(shell pwd)/gobwa/bwa -L$(shell pwd)/jemalloc/lib
export GOPATH=$(shell pwd)

VERSION=0.1

all: arachne gobwa/bwa/libbwa.a gobwa/bwa/bwa

arachne: gobwa/bwa/libbwa.a gobwa/bwa/bwa jemalloc/lib/libjemalloc_pic.a
	@echo "Building arachne"
	mkdir -p bin/
	go build -ldflags "-X arachne/aligner.VERSION=$(VERSION)" -o bin/$@
	cp gobwa/bwa/bwa bin/
	chmod +x bin/arachne

gobwa/bwa/libbwa.a gobwa/bwa/bwa &:
	@echo "Building BWA"
	$(MAKE) -j 4 -C gobwa/bwa libbwa.a bwa

jemalloc/Makefile:
	cd jemalloc && ./autogen.sh && \
	./configure --disable-shared --enable-static

jemalloc/lib/libjemalloc_pic.a: jemalloc/Makefile
	$(MAKE) -j 4 -C jemalloc build_lib_static

clean:
	@echo "Cleaning Build"
	rm -Rf bin/
	$(MAKE) -j 4 -C gobwa/bwa clean
	$(MAKE) -j 4 -C jemalloc distclean
