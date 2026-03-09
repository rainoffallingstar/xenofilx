.PHONY: build test clean install

# Build the xenofilter binary
build:
	go build -o bin/xenofilter ./cmd/xenofilter

# Run tests
test:
	go test -v ./...

# Run tests with coverage
test-coverage:
	go test -coverprofile=coverage.out ./...
	go tool cover -html=coverage.out -o coverage.html

# Install xenofilter to GOPATH/bin
install:
	go install ./cmd/xenofilter

# Clean build artifacts
clean:
	rm -rf bin/
	rm -f coverage.out coverage.html

# Cross-compilation
build-linux:
	GOOS=linux GOARCH=amd64 go build -o bin/xenofilter-linux-amd64 ./cmd/xenofilter

build-mac:
	GOOS=darwin GOARCH=amd64 go build -o bin/xenofilter-darwin-amd64 ./cmd/xenofilter

build-windows:
	GOOS=windows GOARCH=amd64 go build -o bin/xenofilter-windows-amd64.exe ./cmd/xenofilter
