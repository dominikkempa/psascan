#include <cstdio>
#include <cstdlib>
#include <string>

#include "utils.h"
#include "async_stream_reader.h"
#include "async_stream_writer.h"


int main(int argc, char **argv) {
  if (argc != 3) {
    fprintf(stderr, "Usage: %s SIZE OUTFILE\n"
        "Generate Fibonacci word of SIZE MiB and write to OUTFILE.\n", argv[0]);
    std::exit(EXIT_FAILURE);
  }

  std::string output_filename = std::string(argv[2]);
  std::size_t length = (std::size_t)std::atol(argv[1]) << 20;
  fprintf(stderr, "Length = %ld (%.2LfMiB)\n", length, 1.L * length / (1L << 20));
  fprintf(stderr, "Output filename = %s\n", output_filename.c_str());

  async_stream_writer<unsigned char> *initial_writer = new async_stream_writer<unsigned char>(output_filename);
  initial_writer->write('a');
  initial_writer->write('b');
  delete initial_writer;

  std::size_t prev_length = 1;
  std::size_t cur_length = 2;

  std::string temp_filename = output_filename + ".temp";
  while (cur_length < length) {
    async_stream_reader<unsigned char> *reader = new async_stream_reader<unsigned char>(output_filename);
    async_stream_writer<unsigned char> *writer = new async_stream_writer<unsigned char>(temp_filename, "w");

    std::size_t tocopy = std::min(length - cur_length, prev_length);
    for (std::size_t j = 0; j < tocopy; ++j)
      writer->write(reader->read());
    delete writer;
    delete reader;

    async_stream_reader<unsigned char> *reader2 = new async_stream_reader<unsigned char>(temp_filename);
    async_stream_writer<unsigned char> *writer2 = new async_stream_writer<unsigned char>(output_filename, "a");

    while (!reader2->empty())
      writer2->write(reader2->read());

    delete reader2;
    delete writer2;

    prev_length = cur_length;
    cur_length += tocopy;
  }

  if (utils::file_exists(temp_filename))
    utils::file_delete(temp_filename);

  fprintf(stderr, "Finished.\n");
}

