#ifndef CORE_H
#define CORE_H

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <ctype.h>
#include <zlib.h>
#include "KmerIndex.h"
#include "kseq.h"
#include "Threadpool.h"


class Pipeline_worker {
public:
    Pipeline_worker(int dpu_n, KmerIndex *KI, std::vector<std::string> &seqs) : dpu_n(dpu_n), KI(KI), seqs(seqs){

    }

    KmerIndex *KI;
    std::vector<std::string> seqs;
    int dpu_n;

    void map();
};

class Core {
public:
    Core(int pw_n, int dpu_n, KmerIndex *KI, std::string &readFile) : pw_n(pw_n), dpu_n(dpu_n), KI(KI), readFile(readFile) {
        read_rna_file(readFile);

        pool = new ThreadPool(pw_n);
        pw = new Pipeline_worker*[pw_n];
        // calling constructor
        int start = 0;
        int end = -1;
        int n = seqs.size() / pw_n;
        int r = seqs.size() % pw_n;
        for (int i = 0; i < pw_n; i++) {
            // allocate RNA reads to pipeline workers
            start = end + 1;
            end = (i < r) ? (start + n + 1) : (start + n);
            std::vector<std::string>::const_iterator first = seqs.begin() + start;
            std::vector<std::string>::const_iterator last = seqs.begin() + end;
            std::vector<std::string> seq_tmp(first, last);
            pw[i] = new Pipeline_worker(dpu_n, KI, seq_tmp);
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
    std::vector<std::string> seqs;
    ThreadPool *pool;
    Pipeline_worker **pw;

    void read_rna_file(std::string &readFile);
    void allocate_pipeline_worker();
};

#endif 