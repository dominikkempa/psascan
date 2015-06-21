SHELL = /bin/sh

CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -O3 -funroll-loops
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -g2

all: test

test:
	$(CC) $(CFLAGS) -o test main.cpp utils.cpp

clean:
	/bin/rm -f *.o

nuclear:
	/bin/rm -f *.o test
