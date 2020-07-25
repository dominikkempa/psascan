SHELL = /bin/sh

CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -funroll-loops -DNDEBUG -O3 -std=c++0x
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -g2 -std=c++0x

all: rank_test

rank_test:
	$(CC) $(CFLAGS) -o rank_test test.cpp -fopenmp

clean:
	/bin/rm -f *.o

test: rank_test
	@./rank_test

nuclear:
	/bin/rm -f rank_test *.o
