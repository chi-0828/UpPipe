#include "Core.h"

#ifndef KSEQ_INIT_READY
#define KSEQ_INIT_READY
KSEQ_INIT(gzFile, gzread)
#endif

// #ifdef DPU_HEAD
// #define DPU_HEAD 1
// #else
#include <dpu>
using namespace dpu;
// #endif

// pipeline worker 
void Pipeline_worker::map() {
    std::cerr << "[DPU] allocating " << dpu_n << " DPUs" << std::endl;
    auto DPUs = DpuSet::allocate(dpu_n);
}

// class core
void Core::read_rna_file(std::string &readFile) {
    // read input
    std::cerr << "[load] loading RNA read file " << readFile << std::endl;
    // read fasta file  
    gzFile fp = 0;
    kseq_t *seq;
    int l = 0;
    fp = gzopen(readFile.c_str(), "r");
    seq = kseq_init(fp);
    while (true) {
        l = kseq_read(seq);
        if (l <= 0) {
        break;
        }
        seqs.emplace_back(seq->seq.s);
        std::string& str = *seqs.rbegin();
        auto n = str.size();
        for (auto i = 0; i < n; i++) {
            char c = str[i];
            c = ::toupper(c);
            if (c=='U') {
                str[i] = 'T';
            } 
        }
    }
    gzclose(fp);
    fp=0;
} 

void Core::allocate_pipeline_worker() {
    for(int i = 0; i < pw_n; ++i) {
        pool->enqueue(&Pipeline_worker::map, pw[i]);
    }
} 