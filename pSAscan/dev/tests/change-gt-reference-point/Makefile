SHELL = /bin/sh

CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -funroll-loops -DNDEBUG -O3 -pthread -std=c++0x
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -g2 -pthread -std=c++0x

all: test

test:
	$(CC) $(CFLAGS) -o test test.cpp utils.cpp -fopenmp -ldivsufsort

clean:
	/bin/rm -f *.o

nuclear:
	/bin/rm -f test *.o
