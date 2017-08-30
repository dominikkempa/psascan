SHELL = /bin/sh

CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -funroll-loops -DNDEBUG -O3 -pthread -std=c++0x -march=native
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -g2 -std=c++0x -pthread

all: test_sequential test_parallel

test_sequential: main.cpp
	$(CC) $(CFLAGS) -o test_sequential main.cpp utils.cpp -ldivsufsort -ldivsufsort64

test_parallel: main.cpp
	$(CC) $(CFLAGS) -o test_parallel main.cpp utils.cpp -ldivsufsort -ldivsufsort64 -fopenmp

clean:
	/bin/rm -f test_sequential test_parallel *.o
	
nuclear:
	/bin/rm -f test_sequential test_parallel *.o tempfile*
