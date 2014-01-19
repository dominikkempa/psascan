SHELL = /bin/sh

CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -static -funroll-loops -DNDEBUG -O3
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -static -g2

all: fast_rank_test

fast_rank_test: fast_rank_test.cpp
	$(CC) $(CFLAGS) -o fast_rank_test fast_rank_test.cpp

clean:
	/bin/rm -f *.o

test: fast_rank_test
	./fast_rank_test

nuclear:
	/bin/rm -f fast_rank_test *.o
