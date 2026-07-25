.PHONY: build test test-coverage clean install build-linux build-mac build-windows

BINARY := xenofilx
COMMAND := ./cmd/xenofilx
VERSION ?= 0.1.0
VERSION_PACKAGE := github.com/rainoffallingstar/xenofilx/pkg/cli
LDFLAGS := -X $(VERSION_PACKAGE).Version=$(VERSION)

# Build the xenofilx binary
build:
	mkdir -p bin
	go build -ldflags "$(LDFLAGS)" -o bin/$(BINARY) $(COMMAND)

# Run tests
test:
	go test -v ./...

# Run tests with coverage
test-coverage:
	go test -coverprofile=coverage.out ./...
	go tool cover -html=coverage.out -o coverage.html

# Install xenofilx to GOPATH/bin
install:
	go install -ldflags "$(LDFLAGS)" $(COMMAND)

# Clean build artifacts
clean:
	rm -rf bin/
	rm -f coverage.out coverage.html

# Cross-compilation
build-linux:
	mkdir -p bin
	GOOS=linux GOARCH=amd64 go build -ldflags "$(LDFLAGS)" -o bin/$(BINARY)-linux-amd64 $(COMMAND)

build-mac:
	mkdir -p bin
	GOOS=darwin GOARCH=amd64 go build -ldflags "$(LDFLAGS)" -o bin/$(BINARY)-darwin-amd64 $(COMMAND)

build-windows:
	mkdir -p bin
	GOOS=windows GOARCH=amd64 go build -ldflags "$(LDFLAGS)" -o bin/$(BINARY)-windows-amd64.exe $(COMMAND)
