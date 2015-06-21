SHELL = /bin/sh

CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -funroll-loops -DNDEBUG -O3
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -g2 -pthread -std=c++0x

all: gen

gen:
	$(CC) $(CFLAGS) -o gen main.cpp utils.cpp

clean:
	/bin/rm -f *.o

nuclear:
	/bin/rm -f gen *.o
