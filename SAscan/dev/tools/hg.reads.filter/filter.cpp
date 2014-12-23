#include <cstdio>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <string>

int main(int argc, char **argv) {
  if (argc == 1) {
    fprintf(stderr, "usage: %s [FILE]...\n"
                    "Filter FILEs and write to stdout.\n",
                    argv[0]);
    std::exit(EXIT_FAILURE);
  }

  for (long argid = 1; argid < argc; ++argid) {
    std::ifstream fin(argv[argid]);
    std::string fout_filename = std::string(argv[argid]) + ".filtered";
    std::ofstream fout(fout_filename.c_str());
    std::string line;
    long line_cnt = 0;
    while (std::getline(fin, line)) {
      if (line_cnt % 4 == 1)
        fout << line << std::endl;
      ++line_cnt;
    }
  }
}
