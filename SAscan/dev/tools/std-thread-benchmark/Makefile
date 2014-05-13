all: example

example:
	g++ -Wall -Wextra -pedantic -Wshadow -DNDEBUG -O3 -funroll-loops main.cpp -o prog -pthread -std=c++0x

clean:
	rm -rf prog *.o
	
nuclear:
	rm -rf prog *.o
    
