SHELL = /bin/sh

CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -funroll-loops -DNDEBUG -O3
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -g2 -pthread -std=c++0x

all: sa

sa:
	$(CC) $(CFLAGS) -o sa main.cpp utils.cpp -ldivsufsort64

clean:
	/bin/rm -f *.o

nuclear:
	/bin/rm -f sa *.o
