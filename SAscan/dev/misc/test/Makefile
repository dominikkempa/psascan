SHELL = /bin/sh

CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -funroll-loops -DNDEBUG -O3 -pthread -fopenmp -std=c++0x
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -g2 -pthread -fopenmp -std=c++0x

all: test_sascan

test_sascan:
	$(CC) $(CFLAGS) -o test_sascan ../utils.cpp test.cpp -ldivsufsort -ldivsufsort64

clean:
	/bin/rm -f *.o

test: test_sascan
	@./test_sascan

nuclear:
	/bin/rm -f test_sascan *.o
