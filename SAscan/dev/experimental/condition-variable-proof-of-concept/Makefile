SHELL = /bin/sh

CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -funroll-loops -DNDEBUG -O3 -pthread -std=c++0x
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -static -g2

all: program

program: main.cpp
	$(CC) $(CFLAGS) -o program main.cpp

clean:
	/bin/rm -f program *.o
	
nuclear:
	/bin/rm -f program *.o
