#ifndef CORE_H
#define CORE_H

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <mutex>
#include <ctype.h>
#include <zlib.h>
#include <deque>
#include "KmerIndex.h"
#include "kseq.h"
#include "Threadpool.h"
#include "dpu_app/dpu_def.h"

using Read_packet = std::vector<char*>;
using Partial_result =  std::vector<dpu_result>;
Read_packet get();

class Pipeline_worker {
public:
    Pipeline_worker(int dpu_n, KmerIndex *KI) : dpu_n(dpu_n), KI(KI) {
        transfer_to_dpu_buf = new std::deque<Read_packet>(dpu_n, Read_packet(0));
        transfer_from_dpu_buf = new std::vector<Partial_result>(dpu_n);
        partial_result_buf = new std::deque<Partial_result>(dpu_n);
    }

    ~Pipeline_worker() {
        delete transfer_to_dpu_buf;
        delete transfer_from_dpu_buf;
        delete partial_result_buf;
    }

    KmerIndex *KI;
    int dpu_n;
    void map();
    Partial_result compare();

    std::deque<Read_packet> *transfer_to_dpu_buf;
    std::vector<Partial_result> *transfer_from_dpu_buf;
    std::deque<Partial_result> *partial_result_buf;
};

class Core {
public:
    Core(int pw_n, int dpu_n, KmerIndex *KI, std::string &readFile) : pw_n(pw_n), dpu_n(dpu_n), KI(KI), readFile(readFile) {

        pool = new ThreadPool(pw_n);
        pw = new Pipeline_worker*[pw_n];
        // calling constructor
        for (int i = 0; i < pw_n; i++) {
            pw[i] = new Pipeline_worker(dpu_n, KI);
        }

        allocate_pipeline_worker();
    }
    ~Core() {
        delete pool;
    }

    int pw_n; // number of pipeline workers
    int dpu_n; // number of DPUs in a pipeline worker
    KmerIndex *KI;
    std::string readFile;
    ThreadPool *pool;
    Pipeline_worker **pw;

    void read_rna_file(std::string &readFile);
    void allocate_pipeline_worker();
};

#endif 