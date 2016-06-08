#include <cstdio>
#include <cstdlib>
#include <string>

#include "async_stream_reader.h"
#include "async_stream_writer.h"
#include "utils.h"
#include "uint40.h"


template<typename saidx_t>
void truncate(std::string full_sa_filename, std::uint64_t suffix_length) {
  // Compute text length.
  std::uint64_t text_length = utils::file_size(full_sa_filename) / sizeof(saidx_t);

  // Initialize reader of full suffix array.
  typedef async_stream_reader<saidx_t> reader_type;
  reader_type *reader = new reader_type(full_sa_filename);

  // Initialize writer of the output suffix array.
  typedef async_stream_writer<saidx_t> writer_type;
  writer_type *writer = new writer_type();

  // Compute SA subsequence.
  std::uint64_t min_suffix = text_length - suffix_length;
  for (std::uint64_t i = 0; i < text_length; ++i) {
    std::uint64_t sa_i = reader->read();
    if (sa_i >= min_suffix)
      writer->write((saidx_t)(sa_i - min_suffix));
  }

  // Clean up.
  delete writer;
  delete reader;
}

int main(int argc, char **argv) {
  if (argc != 3) {
    fprintf(stderr, "Usage: %s FILE n\n"
        "Compute the suffix array containing the last n suffixes from the\n"
        "suffix array of the full file. The suffix array is read from FILE.\n"
        "The result is written on standard output. The program assumes that\n"
        "suffix array is encoded with 40-bit integers.\n", argv[0]);
    std::exit(EXIT_FAILURE);
  }

  std::uint64_t n = std::atol(argv[2]);
  std::string full_sa_filename = argv[1];
  truncate<uint40>(full_sa_filename, n);
}

