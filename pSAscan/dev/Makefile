SHELL = /bin/sh
CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -funroll-loops -pthread -std=c++0x -DNDEBUG -O3 -march=native
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -pthread -std=c++0x -g2

all: psascan

psascan:
	$(CC) $(CFLAGS) -o psascan src/psascan_src/utils.cpp src/main.cpp -ldivsufsort -ldivsufsort64 -fopenmp

clean:
	/bin/rm -f *.o

nuclear:
	/bin/rm -f psascan *.o
