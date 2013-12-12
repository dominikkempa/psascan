SHELL = /bin/sh

CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -static -funroll-loops -DNDEBUG -O3
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -static -g2

all: fgm

fgm: main.cpp
	$(CC) $(CFLAGS) -o fgm main.cpp

clean:
	/bin/rm -f *.o

test: fgm
	./fgm

nuclear:
	/bin/rm -f fgm *.o gap.* sparseSA.* gt_head new_gt_head gt_tail new_gt_tail
