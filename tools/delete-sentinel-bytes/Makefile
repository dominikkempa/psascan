SHELL = /bin/sh

CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -funroll-loops -DNDEBUG -O3 -march=native -std=c++0x

all: delete_sentinel_bytes

delete_sentinel_bytes:
	$(CC) $(CFLAGS) -o delete_sentinel_bytes ./src/main.cpp

clean:
	/bin/rm -f *.o

nuclear:
	/bin/rm -f delete_sentinel_bytes *.o
