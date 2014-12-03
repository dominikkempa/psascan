SHELL = /bin/sh

CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -funroll-loops -DNDEBUG -O3 -pthread -std=c++0x -msse4.2
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -g2

all: test

test:
	$(CC) $(CFLAGS) -o test main.cpp utils.cpp -fopenmp -ldivsufsort64

clean:
	/bin/rm -f *.o

nuclear:
	/bin/rm -f test *.o
