CC=gcc
CFLAGS=-g -O3 -Wall -Wpedantic -lz -std=gnu11
PROG=fastDBarcode

all:
	mkdir -p ./build
	$(CC) $(CFLAGS) -o ./build/$(PROG) ./src/main.c

clean:
	rm -rvf ./build
