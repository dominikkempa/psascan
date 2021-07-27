pSAscan - Parallel external memory suffix array construction
============================================================


Description
-----------

This repository contains the implementation of pSAscan, a parallel
external-memory suffix array construction algorithm. The basic idea of
the algorithm is to first construct the suffix arrays for blocks of
text and then merge them into the full suffix array.

The key features of pSAscan are:
- The algorithm is able to handle inputs much larger than considered
  in the literature so far. In the experiments, computing the suffix
  array of a 1TiB file with the new algorithm on a machine equipped
  with 120GiB of RAM took a little over a week and required only
  7.2TiB of disk space (including input text and output suffix array),
  whereas on the same machine the previously best algorithm would
  require 3.5 times as much disk space and take about four times
  longer.
- The algorithm is able to fully utilize the available RAM, for
  example, computing the suffix array of a 200GiB file using 3.5GiB of
  RAM takes around 170h. With 120GiB of RAM, the computation time is
  reduced to less than 12h.

For more detailed description of the algorithm and reports of
experimental evaluation, refer to the following paper.

    @inproceedings{kkp15cpm,
      author =    {Juha K{\"{a}}rkk{\"{a}}inen and Dominik Kempa
                   and Simon J. Puglisi},
      title =     {Parallel External Memory Suffix Sorting},
      booktitle = {26th Annual Symposium on Combinatorial Pattern
                   Matching (CPM 2015)},
      pages     = {329--342},
      year      = {2015},
      doi       = {10.1007/978-3-319-19929-0\_28},
    }

The latest version of pSAscan is available from
https://github.com/dominikkempa/psascan.



Compilation and usage
---------------------

1. Download https://github.com/y-256/libdivsufsort/archive/2.0.1.tar.gz
and install. Make sure to compile it to static 64-bit libraries. More
detailed instructions on the (recommended) installation can be found
in the doc/divsufsort-install-guide.txt file of this package.

2. After installing libdivsufsort, pSAscan is compiled by simply
typing `make` in the directory containing this README. This will build
the pSAscan executable called `construct_sa`. For usage instructions,
run the program without any arguments.

### Example

The simplest usage of pSAscan is as follows. Suppose the text is
located in `/data/input.txt`. Then, to compute the suffix array of
`input.txt`, type:

    $ ./construct_sa /data/input.txt

This will write the output suffix array to `/data/input.txt.sa5`. By
default, pSAscan uses 3.5GiB of RAM. The current implementation
encodes the output suffix array using unsigned 40-bit integers. For
further processing of the suffix array, one should use the same or
compatible encoding. The class implementing the unsigned 40-bit
integers is located in the `include/types/uint40.hpp` file.
A more advanced usage of pSAscan is demonstrated below.

    $ ./construct_sa /data/input.txt -m 8gi -o ~/out/sa.out

Explanation:
- The -o flag allows specifying the location and filename of the
  output suffix array. The default location and filename is the same
  as input text, with the appended ".sa5" suffix.
- The -m flag allows specifying the amount of RAM used during the
  computation (in bytes). In this example, the RAM limit is set to 8gi
  = 8 * 2^30 bytes (see the explanation below).

Notes:
- The argument of the -m flag (RAM used during the computation) can be
  specified either explicitly or using common suffixes such as K, M,
  G, T, Ki, Mi, Gi, Ti, which respectively correspond to multipliers:
  10^3, 10^6, 10^9, 10^12, 2^10, 2^20, 2^30, 2^40.  Suffix names are
  not case-sensitive, e.g., Ti = ti, k = K.
- The flags specifying RAM usage, output filename, etc. can be given
  in any order.
- Filenames passed as command-line arguments can be given as relative
  paths, e.g., `../input.txt` and `~/out/sa.out` are valid paths, see
  also example above.



Disk space requirements
-----------------------

To compute the suffix array of an n-byte input text, pSAscan needs
about 7.5n bytes of disk space. This includes the input (n bytes) and
output (5n bytes). In the default mode, pSAscan assumes, that
there is 6.5n bytes of free disk space available in the location used
as the destination for the suffix array. This space is used for
auxiliary files created during the computation and to accommodate the
output.

The above disk space requirement may in some cases prohibit the use of
algorithm, e.g., if there is enough space (5n) on one physical disk to
hold the suffix array, but not enough (6.5n) to run the algorithm. To
still allow the computation in such cases, pSAscan implements the -g
flag. With this flag, one can force pSAscan to use disk space from two
physically different locations (e.g., on two disks). More precisely,
out of 6.5n bytes of disk space used by pSAscan, about n bytes is used
to store the so-called "gap array". By default, the gap array is
stored along with the suffix array. The -g flag allows explicitly
specifying the location of the gap array. This way, it suffices that
there is only 5.5n bytes of disk space in the location specified as
the destination of the suffix array. The remaining n bytes can be
allocated in other location specified with the -g flag.

### Example

Assume the location of input/output files and RAM usage as in the
example from the previous section. To additionally specify the
location of the gap array as `/data3/tmp` run the `construct_sa`
command as:

    $ ./construct_sa /data/input.txt -m 8gi -o /data2/sa.out -g /data3/tmp



RAM requirements
----------------

The algorithm does not have a fixed memory requirements. In principle,
it can run with any amount of RAM (though there is some minimal
per-thread amount necessary in the streaming phase). However, since
the time complexity (without logarithmic factors) of the algorithm is
O(n^2 / M), where M is the amount of RAM used in the computation,
using more RAM decreases the runtime.  Thus, the best performance is
achieved when nearly all unused RAM available in the system (as shown
by the Linux `free` command) is used for the computation. Leaving
about 5% (but not more than 2GiB) of RAM free is advised to prevent
thrashing.

### Example

On a machine with 12 physical cores and Hyper-Threading (and thus
capable of simultaneously running 24 threads) it takes about a week to
compute a suffix array of a 200GiB file using 3.5GiB of RAM. Using
120GiB of RAM reduces the time to less than 12 hours.



Troubleshooting
---------------

1. I am getting an error about the exceeded number of opened files.

Solution: The error is caused by the operating system imposing a limit
on the maximum number of files opened by a program. The limit can be
increased with the `ulimit -n newlimit` command. However, in Linux the
limit cannot be increased beyond the so-called "hard limit", which is
usually only few times larger. Furthermore, this is a temporary
solution that needs to repeated every time a new session is
started. To increase the limits permanently, edit (as a root) the file
`/etc/security/limits.conf` and add the following lines at the end
(including the asterisks):

    * soft nofile 128000
    * hard nofile 128000

This increases the limit to 128000 (use larger values if necessary).
The new limits apply (check with `ulimit -n`) after starting new
session.

2. Program stops without any error message.

Solution: Most likely the problem occurred during internal-memory
sorting.  Re-running the program with the -v flag should show the
error message.



Limitations
-----------

1. The maximum size of input text is 1TiB (2^40 bytes).
2. The current implementation supports only inputs over byte alphabet.
3. Only texts not containing bytes with value 255 are handled
   correctly.  The bytes with value 255 can be removed from the input
   using the tool located in the directory tools/delete-sentinel-bytes/
   of this package.
4. The current internal-memory suffix sorting algorithm used
   internally in pSAscan works only if the input text is split into
   segments of size at most 2GiB each. Therefore, pSAscan will fail,
   if the memory budget X for the computation (specified with the -m
   flag) satisfies X / p > 10 * 2^31, where p is the number of threads
   used during the computation. On most systems, this is not a severe
   limitation, e.g., for a regular 4-core machine supporting
   Hyper-Threading (and thus capable of simultaneously running 8
   threads), pSAscan can utilize up to 160GiB of RAM.

The above limitations (except possibly 2) are not inherent to the
algorithm but rather the current implementation. Future releases will
most likely overcome these limitations.



Third-party code
----------------

This implementation makes use of some third-party code:
- The internal suffix-sorting routine is divsufsort 2.0.1.
  See: https://github.com/y-256/libdivsufsort



Terms of use
------------

pSAscan is released under the MIT/X11 license. See the file LICENCE
for more details. If you use this code, please cite the paper
mentioned above.



Authors
-------

pSAscan was implemented by:
- [Dominik Kempa](https://scholar.google.com/citations?user=r0Kn9IUAAAAJ)
- [Juha Karkkainen](https://scholar.google.com/citations?user=oZepo1cAAAAJ)
