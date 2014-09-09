SHELL = /bin/sh

CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -funroll-loops -DNDEBUG -O3 -pthread -std=c++0x
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -g2 -pthread -std=c++0x

all: test_file test_random bwt_test_file bwt_test_random

test_file:
	$(CC) $(CFLAGS) -o test_file test_file.cpp utils.cpp -fopenmp -ldivsufsort -ldivsufsort64

bwt_test_file:
	$(CC) $(CFLAGS) -o bwt_test_file bwt_test_file.cpp utils.cpp -fopenmp -ldivsufsort -ldivsufsort64

test_random:
	$(CC) $(CFLAGS) -o test_random test_random.cpp utils.cpp -fopenmp -ldivsufsort -ldivsufsort64

bwt_test_random:
	$(CC) $(CFLAGS) -o bwt_test_random bwt_test_random.cpp utils.cpp -fopenmp -ldivsufsort -ldivsufsort64

clean:
	/bin/rm -f *.o

nuclear:
	/bin/rm -f test_file test_random bwt_test_file bwt_test_random *.o
