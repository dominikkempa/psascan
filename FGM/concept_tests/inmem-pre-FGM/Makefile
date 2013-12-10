SHELL = /bin/sh

CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -static -funroll-loops -DNDEBUG -O3
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -static -g2

all: inmem_fgm

inmem_fgm: inmem_fgm.cpp
	$(CC) $(CFLAGS) -o inmem_fgm inmem_fgm.cpp

clean:
	/bin/rm -f *.o

test: inmem_fgm
	./inmem_fgm

nuclear:
	/bin/rm -f inmem_fgm *.o
