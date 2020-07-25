SHELL = /bin/sh

CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -funroll-loops -pthread -std=c++0x -DNDEBUG -O3
#CFLAGS = -Wall -Wextra -Wshadow -pedantic -g2 -std=c++0x -pthread

all: test test_file test_file_with_sa

test:
	$(CC) $(CFLAGS) -o test test.cpp utils.cpp -ldivsufsort -fopenmp

test_file:
	$(CC) $(CFLAGS) -o test_file test_file.cpp utils.cpp -ldivsufsort -fopenmp

test_file_with_sa:
	$(CC) $(CFLAGS) -o test_file_with_sa test_file_with_sa.cpp utils.cpp -ldivsufsort -fopenmp

clean:
	/bin/rm -f *.o

nuclear:
	/bin/rm -f test test_file test_file_with_sa *.o
