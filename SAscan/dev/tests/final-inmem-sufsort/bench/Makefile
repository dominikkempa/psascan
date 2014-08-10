SHELL = /bin/sh

CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -funroll-loops -DNDEBUG -O3 -pthread -std=c++0x
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -g2 -pthread -std=c++0x

all: test_file

test_file:
	$(CC) $(CFLAGS) -o test_file test_file.cpp utils.cpp -fopenmp -ldivsufsort

clean:
	/bin/rm -f *.o

nuclear:
	/bin/rm -f test_file *.o
