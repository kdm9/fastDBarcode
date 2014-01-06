CC=gcc
DEBUG_FLAGS=-g -pg
CFLAGS=$(DEBUG_FLAGS) -O3 -Wall -Wpedantic -lz -std=gnu11 -fopenmp
PROG=fastDBarcode

all:
	mkdir -p ./build
	$(CC) $(CFLAGS) -o ./build/$(PROG) ./src/main.c

clean:
	rm -rvf ./build
