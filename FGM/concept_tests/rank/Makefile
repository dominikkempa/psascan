SHELL = /bin/sh

CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -static -funroll-loops -DNDEBUG -O3
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -static -g2

all: general_rank_test

general_rank_test: general_rank_test.cpp
	$(CC) $(CFLAGS) -o general_rank_test general_rank_test.cpp

clean:
	/bin/rm -f *.o

test: general_rank_test
	./general_rank_test

nuclear:
	/bin/rm -f general_rank_test *.o
