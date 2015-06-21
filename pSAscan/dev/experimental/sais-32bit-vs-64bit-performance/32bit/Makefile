SHELL = /bin/sh

CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -static -funroll-loops -DNDEBUG -O3
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -static -g2

all: sufsort

sufsort: main.cpp utils.cpp
	$(CC) $(CFLAGS) -o sufsort utils.cpp main.cpp

clean:
	/bin/rm -f *.o

nuclear:
	/bin/rm -f sufsort *.o
