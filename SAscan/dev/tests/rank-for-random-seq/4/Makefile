SHELL = /bin/sh

CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -funroll-loops -DNDEBUG -O3 -std=c++0x
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -g2 -std=c++0x

all: rank_test test_file

rank_test:
	$(CC) $(CFLAGS) -o rank_test test.cpp -fopenmp

test_file:
	$(CC) $(CFLAGS) -o test_file test_file.cpp -fopenmp

clean:
	/bin/rm -f *.o

test: rank_test
	@./rank_test

nuclear:
	/bin/rm -f rank_test test_file *.o
