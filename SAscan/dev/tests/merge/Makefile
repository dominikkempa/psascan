SHELL = /bin/sh

CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -static -funroll-loops -DNDEBUG -O3
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -static -g2

all: fast_merge_test fast_merge_test_2

fast_merge_test: fast_merge_test.cpp
	$(CC) $(CFLAGS) -o fast_merge_test fast_merge_test.cpp

fast_merge_test_2: fast_merge_test_2.cpp
	$(CC) $(CFLAGS) -o fast_merge_test_2 fast_merge_test_2.cpp

clean:
	/bin/rm -f *.o

test: fast_merge_test fast_merge_test_2
	@./fast_merge_test
	@./fast_merge_test_2


nuclear:
	/bin/rm -f fast_merge_test fast_merge_test_2 *.o
