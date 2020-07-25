SHELL = /bin/sh

CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -static -funroll-loops -DNDEBUG -O3
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -static -g2

all: test_fast_rank

test_fast_rank: test_fast_rank.cpp
	$(CC) $(CFLAGS) -o test_fast_rank test_fast_rank.cpp

clean:
	/bin/rm -f *.o

test: test_fast_rank
	@./test_fast_rank

nuclear:
	/bin/rm -f test_fast_rank *.o
