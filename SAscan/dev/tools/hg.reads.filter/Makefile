SHELL = /bin/sh

CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -static -funroll-loops -DNDEBUG -O3

all: filter

filter: filter.cpp
	$(CC) $(CFLAGS) -o filter filter.cpp

clean:
	/bin/rm -f filter *.o
	
nuclear:
	/bin/rm -f filter *.o
