SHELL = /bin/sh

CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -funroll-loops -DNDEBUG -O3 -pthread -std=c++0x -march=native
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -static -g2

all: test

test: main.cpp
	$(CC) $(CFLAGS) -o test main.cpp utils.cpp -ldivsufsort -ldivsufsort64

clean:
	/bin/rm -f test *.o
	
nuclear:
	/bin/rm -f test *.o tempfile*
