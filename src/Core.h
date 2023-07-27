#ifndef CORE_H
#define CORE_H

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <thread>
#include <mutex>
#include <ctype.h>
#include <zlib.h>
#include <deque>
#include "KmerIndex.h"
#include "kseq.h"
#include "dpu_app/dpu_def.h"

using Read_packet = std::vector<char>;
using Partial_result =  std::vector<dpu_result>;
Read_packet get();

class Pipeline_worker {
public:
    
    ~Pipeline_worker() {
        delete transfer_to_dpu_buf;
        delete transfer_from_dpu_buf;
        delete partial_result_buf;
    }

    KmerIndex *KI;
    int dpu_n;
    void operator()(int n, KmerIndex *ki);
    void map();
    Partial_result compare();

    std::deque<Read_packet> *transfer_to_dpu_buf;
    std::vector<Partial_result> *transfer_from_dpu_buf;
    std::deque<Partial_result> *partial_result_buf;
};

class Core {
public:
    Core(int pw_n, int dpu_n, KmerIndex *KI, std::string &readFile) : pw_n(pw_n), dpu_n(dpu_n), KI(KI), readFile(readFile) {
        std::cerr << "[UpPipe] start alignment" << std::endl;
        allocate_pipeline_worker();
    }
    ~Core() {
        std::cerr << "[UpPipe] end alignment" << std::endl;
    }

    int pw_n; // number of pipeline workers
    int dpu_n; // number of DPUs in a pipeline worker
    KmerIndex *KI;
    std::string readFile;

    void read_rna_file(std::string &readFile);
    void allocate_pipeline_worker();
};

#endif 