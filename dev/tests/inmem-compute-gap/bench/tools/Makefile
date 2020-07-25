SHELL = /bin/sh

CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -funroll-loops -pthread -std=c++0x -DNDEBUG -O3
#CFLAGS = -Wall -Wextra -Wshadow -pedantic -g2 -std=c++0x -pthread

all: compute_sa_to_file

compute_sa_to_file:
	$(CC) $(CFLAGS) -o compute_sa_to_file main.cpp utils.cpp -ldivsufsort -fopenmp

clean:
	/bin/rm -f *.o

nuclear:
	/bin/rm -f test compute_sa_to_file *.o
