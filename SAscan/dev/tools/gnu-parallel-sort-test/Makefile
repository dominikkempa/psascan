SHELL = /bin/sh

CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -funroll-loops -static -fopenmp -DNDEBUG -O3
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -pthread -std=c++0x -g2

all: test

test:
	$(CC) $(CFLAGS) -o test main.cpp

clean:
	/bin/rm -f *.o

nuclear:
	/bin/rm -f test
