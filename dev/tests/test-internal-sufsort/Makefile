SHELL = /bin/sh

CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -funroll-loops -pthread -std=c++0x -DNDEBUG -O3
#CFLAGS = -Wall -Wextra -Wshadow -pedantic -g2 -std=c++0x -pthread

all: test bench

test:
	$(CC) $(CFLAGS) -o test test_internal_sufsort.cpp -ldivsufsort -ldivsufsort64 -fopenmp

bench:
	$(CC) $(CFLAGS) -o bench bench_internal_sufsort.cpp -ldivsufsort -ldivsufsort64 -fopenmp

clean:
	/bin/rm -f *.o

nuclear:
	/bin/rm -f bench test *.o
