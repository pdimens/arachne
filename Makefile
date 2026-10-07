export CGO_LDFLAGS = -L$(shell pwd)/gominibwa/minibwa

VERSION=0.2.1

all: arachne gominibwa/minibwa/libminibwa.a gominibwa/minibwa/minibwa

arachne: gominibwa/minibwa/libminibwa.a gominibwa/minibwa/minibwa
	@echo "Building arachne"
	mkdir -p bin/
	go build -ldflags "-X arachne/aligner.VERSION=$(VERSION) -s -w" -o bin/$@
	cp gominibwa/minibwa/minibwa bin/
	chmod +x bin/arachne

gominibwa/minibwa/libminibwa.a gominibwa/minibwa/minibwa &:
	@echo "Building minibwa"
	$(MAKE) -j 4 -C gominibwa/minibwa libminibwa.a minibwa

clean:
	@echo "Cleaning Build"
	rm -Rf bin/
	$(MAKE) -j 4 -C gominibwa/minibwa clean
