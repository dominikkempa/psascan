all: example

example:
	g++ -Wall -Wextra -pedantic -O3 -DNDEBUG -funroll-loops -static -fopenmp main.cpp -o prog

clean:
	rm -rf prog *.o
	
nuclear:
	rm -rf prog *.o
    
