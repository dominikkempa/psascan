SHELL = /bin/sh

CC = g++
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -funroll-loops -DNDEBUG -O3 -pthread -std=c++0x -march=native
CFLAGS = -Wall -Wextra -pedantic -Wshadow -g2 -pthread -std=c++0x
AUX_PAR_FLAGS = -fopenmp

all: test_parallel test_sequential

test_parallel:
	$(CC) $(CFLAGS) -o test_parallel main.cpp utils.cpp $(AUX_PAR_FLAGS)

test_sequential:
	$(CC) $(CFLAGS) -o test_sequential main.cpp utils.cpp

clean:
	/bin/rm -f *.o

nuclear:
	/bin/rm -f test_parallel test_sequential *.o
