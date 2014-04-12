SHELL = /bin/sh

CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -static -funroll-loops -DNDEBUG -O3
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -static -g2

all: sascan

sascan: main.cpp utils.cpp
	$(CC) $(CFLAGS) -o sascan utils.cpp main.cpp -ldivsufsort -ldivsufsort64

clean:
	/bin/rm -f *.o

test_sascan:
	@cd test; make test

nuclear:
	/bin/rm -f sascan *.o gap.* sparseSA.* gt_head new_gt_head gt_tail new_gt_tail
