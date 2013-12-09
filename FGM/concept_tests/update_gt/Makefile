SHELL = /bin/sh

CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -static -funroll-loops -DNDEBUG -O3
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -static -g2

all: updating_gt_test

updating_gt_test: updating_gt_test.cpp
	$(CC) $(CFLAGS) -o updating_gt_test updating_gt_test.cpp

clean:
	/bin/rm -f *.o

test: updating_gt_test
	./updating_gt_test

nuclear:
	/bin/rm -f updating_gt_test *.o
