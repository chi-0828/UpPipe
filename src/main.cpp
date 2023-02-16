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
#include "KmerIndex.h"
#include "common.h"
#include <cstdio>
#include <zlib.h>

using namespace std;


int main(int argc, char** argv) {
    ProgramOptions opt;
    if(argv[1] == "aligment") {
        ParseOptionsAligment(argc-1,argv+1,opt);
        //KmerIndex *KI = new KmerIndex(opt);
    }
    else {
        ParseOptionsBuild(argc-1,argv+1,opt);
        KmerIndex *KI = new KmerIndex(opt);
        KI->Build(opt);
    }
}