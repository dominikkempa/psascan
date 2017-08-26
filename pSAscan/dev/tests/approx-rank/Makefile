SHELL = /bin/sh

CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -funroll-loops -DNDEBUG -O3 -pthread -std=c++0x -march=native -fopenmp
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -g2 -pthread -std=c++0x -fopenmp

all: test

test:
	$(CC) $(CFLAGS) -o test main.cpp utils.cpp

clean:
	/bin/rm -f *.o

nuclear:
	/bin/rm -f test *.o
