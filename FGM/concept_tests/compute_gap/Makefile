SHELL = /bin/sh

CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -static -funroll-loops -DNDEBUG -O3
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -static -g2

all: gap_array_test

gap_array_test: gap_array_test.cpp
	$(CC) $(CFLAGS) -o gap_array_test gap_array_test.cpp

clean:
	/bin/rm -f *.o

test: gap_array_test
	./gap_array_test

nuclear:
	/bin/rm -f gap_array_test *.o
