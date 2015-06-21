SHELL = /bin/sh

CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -O3 -funroll-loops
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -g2

all: gen del

del: main_del.cpp
	$(CC) $(CFLAGS) -o del main_del.cpp

gen: main_gen.cpp
	$(CC) $(CFLAGS) -o gen main_gen.cpp

clean:
	/bin/rm -f *.o

nuclear:
	/bin/rm -f *.o del gen
