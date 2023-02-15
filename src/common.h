#ifndef KALLISTO_COMMON_H
#define KALLISTO_COMMON_H

#define KALLISTO_VERSION "0.48.0"

#include <string>
#include <vector>
#include <iostream>
#include <getopt.h>

struct ProgramOptions {
  int dpu_n;
  int worker_n;
  std::string index;
  int k;
  std::string output;
  std::vector<std::string> transfasta;
  bool single_end;
  std::string chromFile;
  std::string ecFile;
  std::string transcriptsFile;
  std::string genemap;

ProgramOptions() :
  k(31),
  dpu_n(60),
  worker_n(32),
  single_end(true)
  {}
};

std::string pretty_num(size_t num);
std::string pretty_num(int64_t num);
std::string pretty_num(int num);

void ParseOptionsAligment(int argc, char **argv, ProgramOptions& opt);
void ParseOptionsBuild(int argc, char **argv, ProgramOptions& opt);


#endif // KALLISTO_COMMON_H
