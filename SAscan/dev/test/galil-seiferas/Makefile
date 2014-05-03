SHELL = /bin/sh

CC = g++
CPPFLAGS = -Wall -Wextra -Wshadow -pedantic -funroll-loops -O3 -DNDEBUG

all: strrank_test GS_sets_test smaller_suffixes_inmem_test

strrank_test:
	$(CC) $(CPPFLAGS) -o strrank_test strrank_test.cpp srank.cpp divsufsort.c

smaller_suffixes_inmem_test:
	$(CC) $(CPPFLAGS) -o smaller_suffixes_inmem_test smaller_suffixes_inmem_test.cpp srank.cpp divsufsort.c

GS_sets_test:
	$(CC) $(CPPFLAGS) -o GS_sets_test GS_sets_test.cpp srank.cpp divsufsort.c

clean:
	rm -f strrank_test GS_sets_test smaller_suffixes_inmem_test

nuclear:
	rm -f strrank_test GS_sets_test smaller_suffixes_inmem_test *.o
