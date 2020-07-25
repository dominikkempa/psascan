SHELL = /bin/sh

CC = g++
CPPFLAGS = -Wall -Wextra -Wshadow -pedantic -funroll-loops -O3 -DNDEBUG

all: test

test:
	$(CC) $(CPPFLAGS) -o test main.cpp utils.cpp smaller_suffixes.cpp

clean:
	rm -f test

nuclear:
	rm -f test *.o text text.gt_tail
