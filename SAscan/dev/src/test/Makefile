SHELL = /bin/sh

CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -funroll-loops -DNDEBUG -O3 -pthread -std=c++0x
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -g2 -pthread -std=c++0x

all: test_sascan

test_sascan: test.cpp ../utils.cpp
	$(CC) $(CFLAGS) -o test_sascan ../utils.cpp ../smaller_suffixes.cpp test.cpp -ldivsufsort -ldivsufsort64

clean:
	/bin/rm -f *.o

test: test_sascan
	@./test_sascan

nuclear:
	/bin/rm -f test_sascan *.o gap.* sparseSA.* gt_head new_gt_head gt_tail new_gt_tail
