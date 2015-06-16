SHELL = /bin/sh

CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -funroll-loops -DNDEBUG -O3 -pthread -std=c++0x -march=native
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -g2 -std=c++0x -pthread

all: test

test: main.cpp
	$(CC) $(CFLAGS) -o test main.cpp utils.cpp -ldivsufsort -ldivsufsort64 -fopenmp

clean:
	/bin/rm -f test *.o
	
nuclear:
	/bin/rm -f test *.o tempfile*
