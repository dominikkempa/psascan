SHELL = /bin/sh

CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -static -funroll-loops -DNDEBUG -O3
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -static -g2

all: fgm_test

fgm_test: main.cpp
	$(CC) $(CFLAGS) -o fgm_test main.cpp

clean:
	/bin/rm -f *.o

test: fgm_test
	./fgm_test

nuclear:
	/bin/rm -f fgm_test *.o gap.* sparseSA.* gt_head new_gt_head gt_tail new_gt_tail
