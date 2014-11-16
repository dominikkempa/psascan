SHELL = /bin/sh

CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -funroll-loops -DNDEBUG -O3 -std=c++0x
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -static -g2

all: test_distrib_file

test_distrib_file: main.cpp
	$(CC) $(CFLAGS) -o test_distrib_file main.cpp utils.cpp -fopenmp

clean:
	/bin/rm -f *.o

test: test_distrib_file
	@./test_distrib_file


nuclear:
	/bin/rm -f test_distrib_file *.o temp*
