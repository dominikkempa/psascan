SHELL = /bin/sh

CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -static -funroll-loops -DNDEBUG -O3
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -static -g2

all: sa_correcting_test sa_correcting_test_2 sa_correcting_test_3

sa_correcting_test: sa_correcting_test.cpp
	$(CC) $(CFLAGS) -o sa_correcting_test sa_correcting_test.cpp

sa_correcting_test_2: sa_correcting_test_2.cpp
	$(CC) $(CFLAGS) -o sa_correcting_test_2 sa_correcting_test_2.cpp

sa_correcting_test_3: sa_correcting_test_3.cpp
	$(CC) $(CFLAGS) -o sa_correcting_test_3 sa_correcting_test_3.cpp

clean:
	/bin/rm -f *.o

test: sa_correcting_test sa_correcting_test_2 sa_correcting_test_3
	./sa_correcting_test
	./sa_correcting_test_2
	./sa_correcting_test_3


nuclear:
	/bin/rm -f sa_correcting_test sa_correcting_test_2 sa_correcting_test_3 *.o
