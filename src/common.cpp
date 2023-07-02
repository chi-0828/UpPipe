#include "common.h"

std::string pretty_num(int num) {
  return pretty_num(static_cast<size_t>(num));
}

std::string pretty_num(int64_t num) {
  if (num < 0) {
    return "-" + pretty_num(static_cast<size_t>(num));
  } else {
    return pretty_num(static_cast<size_t>(num));
  }
}

std::string pretty_num(size_t num) {
  auto s = std::to_string(num);
  auto ret = std::string("");

  if (s.size() <= 3) {
    return s;
  }

  int remainder = s.size() % 3;
  if (remainder == 0) {
    remainder = 3;
  }

  size_t start_pos = 0;
  while (start_pos + remainder < s.size() - 1) {
    ret += s.substr(start_pos, remainder) + ",";
    start_pos += remainder;
    remainder = 3;
  }

  ret += s.substr(start_pos, 3);

  return ret;
}

void ParseOptionsAligment(int argc, char **argv, ProgramOptions& opt) {
  const char *opt_string = "i:d:r:o:f:";
  static struct option long_options[] = {
    // long args
    {"index", 1, NULL, 'i'},
    {"dpu", 1, NULL, 'd'},
    {"worker", 1, NULL, 'r'},
    {"output", 1, NULL, 'o'},
    {"read_file", 1, NULL, 'f'},
    { NULL, 0, NULL, 0}
  };
  int c;
  int option_index = 0;
  while (true) {
    c = getopt_long(argc, argv, opt_string, long_options, NULL);
    if (c == -1) {
        break;
    }

    switch (c) {
    case 0:
        break;
    case 'i': {
        opt.index = optarg;
        break;
    }
    case 'd': {
        opt.dpu_n = atoi(optarg);
        break;
    }
    case 'r': {
        opt.worker_n = atoi(optarg);
        break;
    }
    case 'o': {
        opt.output = optarg;
        break;
    }
    case 'f': {
        opt.readFile = optarg;
        break;
    }
    default: 
        break;
    }
  }
}

void ParseOptionsBuild(int argc, char **argv, ProgramOptions& opt) {
  const char *opt_string = "i:d:k:f:";
  static struct option long_options[] = {
    // long args
    {"index", 1, NULL, 'i'},
    {"dpu", 1, NULL, 'd'},
    {"kmer", 1, NULL, 'k'},
    {"read_file", 1, NULL, 'f'},
    {NULL, 0, NULL, 0}
  };
  int c;
  int option_index = 0;
  while (true) {
    c = getopt_long(argc, argv, opt_string, long_options, NULL);
    //std::cerr << (char)c << " " << optarg << "\n";
    if (c == -1) {
        break;
    }

    switch (c) {
    case 0:
        break;
    case 'i': {
        opt.index = optarg;
        break;
    }
    case 'd': {
        opt.dpu_n = atoi(optarg);
        break;
    }
    case 'k': {
        opt.k = atoi(optarg);
        break;
    }
    case 'f': {
        opt.transfasta = optarg;
        break;
    }
    default: 
        break;
    }
  }
}