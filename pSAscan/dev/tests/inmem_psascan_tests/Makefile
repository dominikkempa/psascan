BSHELL = /bin/sh

CC = $(CXX)
CFLAGS = -Wall -Wextra -pedantic -Wshadow -funroll-loops -DNDEBUG -O3 -pthread -std=c++0x -march=native
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -g2 -pthread -std=c++0x

all: sa_test_file sa_test_random psa_test_file psa_test_random bwt_test_file bwt_test_random pbwt_test_file pbwt_test_random gt_test_file gt_test_random pgt_test_file pgt_test_random multi_test_file

sa_test_file:
	$(CC) $(CFLAGS) -o sa_test_file sa_test_file.cpp utils.cpp ../../src/psascan_src/utils.cpp -fopenmp -ldivsufsort -ldivsufsort64

psa_test_random:
	$(CC) $(CFLAGS) -o psa_test_random psa_test_random.cpp utils.cpp ../../src/psascan_src/utils.cpp -fopenmp -ldivsufsort -ldivsufsort64

pbwt_test_random:
	$(CC) $(CFLAGS) -o pbwt_test_random pbwt_test_random.cpp utils.cpp ../../src/psascan_src/utils.cpp -fopenmp -ldivsufsort -ldivsufsort64

pgt_test_file:
	$(CC) $(CFLAGS) -o pgt_test_file pgt_test_file.cpp utils.cpp ../../src/psascan_src/utils.cpp -fopenmp -ldivsufsort -ldivsufsort64

pgt_test_random:
	$(CC) $(CFLAGS) -o pgt_test_random pgt_test_random.cpp utils.cpp ../../src/psascan_src/utils.cpp -fopenmp -ldivsufsort -ldivsufsort64

pbwt_test_file:
	$(CC) $(CFLAGS) -o pbwt_test_file pbwt_test_file.cpp utils.cpp ../../src/psascan_src/utils.cpp -fopenmp -ldivsufsort -ldivsufsort64

psa_test_file:
	$(CC) $(CFLAGS) -o psa_test_file psa_test_file.cpp utils.cpp ../../src/psascan_src/utils.cpp -fopenmp -ldivsufsort -ldivsufsort64

multi_test_file:
	$(CC) $(CFLAGS) -o multi_test_file multi_test_file.cpp utils.cpp ../../src/psascan_src/utils.cpp -fopenmp -ldivsufsort -ldivsufsort64

bwt_test_file:
	$(CC) $(CFLAGS) -o bwt_test_file bwt_test_file.cpp utils.cpp ../../src/psascan_src/utils.cpp -fopenmp -ldivsufsort -ldivsufsort64

gt_test_file:
	$(CC) $(CFLAGS) -o gt_test_file gt_test_file.cpp utils.cpp ../../src/psascan_src/utils.cpp -fopenmp -ldivsufsort -ldivsufsort64

sa_test_random:
	$(CC) $(CFLAGS) -o sa_test_random sa_test_random.cpp utils.cpp ../../src/psascan_src/utils.cpp -fopenmp -ldivsufsort -ldivsufsort64

bwt_test_random:
	$(CC) $(CFLAGS) -o bwt_test_random bwt_test_random.cpp utils.cpp ../../src/psascan_src/utils.cpp -fopenmp -ldivsufsort -ldivsufsort64

gt_test_random:
	$(CC) $(CFLAGS) -o gt_test_random gt_test_random.cpp utils.cpp ../../src/psascan_src/utils.cpp -fopenmp -ldivsufsort -ldivsufsort64

clean:
	/bin/rm -f *.o

nuclear:
	/bin/rm -f sa_test_file sa_test_random bwt_test_file bwt_test_random gt_test_file gt_test_random multi_test_file psa_test_random psa_test_file pbwt_test_file pgt_test_file pbwt_test_random pgt_test_random *.o gt_begin_reversed_bv* supertext.txt
