SHELL = /bin/sh

CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -static -funroll-loops -DNDEBUG -O3
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -static -g2

all: fgm test_fgm

fgm: main.cpp fgm.cpp utils.cpp
	$(CC) $(CFLAGS) -o fgm utils.cpp main.cpp fgm.cpp -ldivsufsort -ldivsufsort64

test_fgm: test.cpp fgm.cpp utils.cpp
	$(CC) $(CFLAGS) -o test_fgm utils.cpp test.cpp fgm.cpp -ldivsufsort -ldivsufsort64

clean:
	/bin/rm -f *.o

test: fgm
	./test_fgm

nuclear:
	/bin/rm -f fgm test_fgm *.o gap.* sparseSA.* gt_head new_gt_head gt_tail new_gt_tail
