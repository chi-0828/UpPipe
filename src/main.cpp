#include <string>
#include <iostream>
#include <random>
#include <sstream>
#include <vector>
#include <sys/stat.h>

#include <thread>
#include <time.h>
#include <algorithm>
#include <limits>
#include <cstdio>
#include <zlib.h>
#include "KmerIndex.h"
#include "common.h"
#include "Core.h"

using namespace std;


int main(int argc, char** argv) {
    ProgramOptions opt;
    if(argc < 2) {
        std::cerr << "Please choose the option: UpPipe {build, alignment} args\n";
        return 0;
    }
    if(strcmp(argv[1], "alignment") == 0) {
        ParseOptionsAligment(argc-1,argv+1,opt);
        KmerIndex *KI = new KmerIndex(opt);
        KI->load(opt);

        Core* core = new Core(opt.worker_n, KI->dpu_n, KI, opt.readFile);
        
        delete core;
        delete KI;
    }
    else if (strcmp(argv[1], "build") == 0) {
        ParseOptionsBuild(argc-1,argv+1,opt);
        KmerIndex *KI = new KmerIndex(opt);
        KI->Build(opt);
        delete KI;
    }
    else {
        std::cerr << "Please choose the option: UpPipe {build, alignment} args\n";
    }

    return 0;
}