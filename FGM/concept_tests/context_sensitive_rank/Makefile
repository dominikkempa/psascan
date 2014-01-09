SHELL = /bin/sh

CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -static -funroll-loops -DNDEBUG -O3
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -static -g2

all: context_rank_test

context_rank_test: context_rank_test.cpp
	$(CC) $(CFLAGS) -o context_rank_test context_rank_test.cpp

clean:
	/bin/rm -f *.o

test: context_rank_test
	./context_rank_test

nuclear:
	/bin/rm -f context_rank_test *.o
