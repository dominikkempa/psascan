all: test

test:
	g++ -Wall -Wextra -pedantic -Wshadow -DNDEBUG -O3 -funroll-loops main.cpp utils.cpp -o test -pthread -std=c++0x

clean:
	rm -rf *.o
	
nuclear:
	rm -rf test *.o

