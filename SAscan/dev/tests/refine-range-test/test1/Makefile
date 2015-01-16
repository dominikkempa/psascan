SHELL = /bin/sh

CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -funroll-loops -DNDEBUG -O3 -pthread -std=c++0x
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -static -g2

all: test

test: main.cpp
	$(CC) $(CFLAGS) -o test main.cpp utils.cpp

clean:
	/bin/rm -f test *.o
	
nuclear:
	/bin/rm -f test *.o
