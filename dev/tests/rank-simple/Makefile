SHELL = /bin/sh

CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -static -funroll-loops -DNDEBUG -O3
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -static -g2

all: simple_rank_test

simple_rank_test:
	$(CC) $(CFLAGS) -o simple_rank_test simple_rank_test.cpp

clean:
	/bin/rm -f *.o

test: simple_rank_test
	@./simple_rank_test

nuclear:
	/bin/rm -f simple_rank_test *.o
