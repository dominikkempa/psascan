SHELL = /bin/sh

CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -static -funroll-loops -DNDEBUG -O3
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -static -g2

all: check_sa

check_sa: check_sa.cpp
	$(CC) $(CFLAGS) -o check_sa check_sa.cpp

clean:
	/bin/rm -f *.o

nuclear:
	/bin/rm -f check_sa *.o
