BSHELL = /bin/sh

CC = $(CXX)
CFLAGS = -Wall -Wextra -pedantic -Wshadow -funroll-loops -DNDEBUG -O2 -pthread -std=c++0x -march=native

all: bench

bench:
	$(CC) $(CFLAGS) -o bench main.cpp utils.cpp ../src/psascan_src/utils.cpp -fopenmp -ldivsufsort -ldivsufsort64

clean:
	/bin/rm -f *.o

nuclear:
	/bin/rm -f *.o bench
