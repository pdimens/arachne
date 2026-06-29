export CGO_LDFLAGS = -L$(shell pwd)/gobwa/bwa -L$(shell pwd)/jemalloc/lib
export GOPATH=$(shell pwd)

VERSION=1.0-dev

all: arachne gobwa/bwa/libbwa.a gobwa/bwa/bwa

arachne: gobwa/bwa/libbwa.a gobwa/bwa/bwa jemalloc/lib/libjemalloc_pic.a
	@echo "Building arachne"
	mkdir -p bin/
	go build -o bin/$@
	cp gobwa/bwa/bwa bin/
	chmod +x bin/arachne

gobwa/bwa/libbwa.a gobwa/bwa/bwa &:
	@echo "Building BWA"
	$(MAKE) -C gobwa/bwa libbwa.a bwa

jemalloc/Makefile:
		cd jemalloc && ./autogen.sh && \
		./configure --disable-shared --enable-static

jemalloc/lib/libjemalloc_pic.a: jemalloc/Makefile
		$(MAKE) -C jemalloc build_lib_static

clean:
	@echo "Cleaning Build"
	rm -Rf bin/
	$(MAKE) -C gobwa/bwa clean
	$(MAKE) -C jemalloc distclean
