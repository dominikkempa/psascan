SHELL = /bin/sh

CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -static -funroll-loops -DNDEBUG -O3
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -static -g2

all: sascan test_sascan

sascan: main.cpp sascan.cpp utils.cpp
	$(CC) $(CFLAGS) -o sascan utils.cpp main.cpp sascan.cpp -ldivsufsort -ldivsufsort64

test_sascan: test.cpp sascan.cpp utils.cpp
	$(CC) $(CFLAGS) -o test_sascan utils.cpp test.cpp sascan.cpp -ldivsufsort -ldivsufsort64

clean:
	/bin/rm -f *.o

test: test_sascan
	./test_sascan

nuclear:
	/bin/rm -f sascan test_sascan *.o gap.* sparseSA.* gt_head new_gt_head gt_tail new_gt_tail
