SHELL = /bin/sh

CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -static -funroll-loops -DNDEBUG -O3
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -static -g2

all: gt_eof_test

gt_eof_test: gt_eof_test.cpp
	$(CC) $(CFLAGS) -o gt_eof_test gt_eof_test.cpp

clean:
	/bin/rm -f *.o

test: gt_eof_test
	@./gt_eof_test

nuclear:
	/bin/rm -f gt_eof_test *.o
