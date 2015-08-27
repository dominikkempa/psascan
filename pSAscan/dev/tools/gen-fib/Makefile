SHELL = /bin/sh

CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -funroll-loops -DNDEBUG -O3 -std=c++0x -pthread
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -g2  -std=c++0x -pthread

all: gen

gen:
	$(CC) $(CFLAGS) -o gen main.cpp utils.cpp

clean:
	/bin/rm -f *.o

nuclear:
	/bin/rm -f gen *.o
