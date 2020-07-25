SHELL = /bin/sh
CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -funroll-loops -std=c++11 -O3 -DNDEBUG -pthread -march=native
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -std=c++11 -g2 -O3 -pthread

all: truncate

truncate:
	$(CC) $(CFLAGS) -o truncate main.cpp utils.cpp

clean:
	/bin/rm -f *.o

nuclear:
	/bin/rm -f truncate *.o
