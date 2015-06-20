SHELL = /bin/sh

CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -funroll-loops -pthread -std=c++0x -DNDEBUG -O3 -march=native
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -pthread -std=c++0x -g2

all: sascan

sascan: main.cpp utils.cpp
	$(CC) $(CFLAGS) -o sascan utils.cpp main.cpp ./inmem_sascan/inmem_sascan_src/utils.cpp -ldivsufsort -ldivsufsort64 -fopenmp

clean:
	/bin/rm -f *.o

test_sascan:
	@cd test; make test

nuclear:
	/bin/rm -f sascan *.o
