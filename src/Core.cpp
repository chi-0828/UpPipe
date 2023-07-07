#include "Core.h"
#include "Dqueue.h"

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

// global RNA read pool
// global RNA read pool
bool end = false;
TDeque RNA_read_pool;
Read_packet get() {
    Read_packet rp;

    while(RNA_read_pool.empty() && !RNA_read_pool.get_end()) ;
    
    
    if(RNA_read_pool.get_end())
        rp = Read_packet(0);
    else
        rp = RNA_read_pool.pop_front();

    return rp;
}

// pipeline worker 
void Pipeline_worker::map() {
    std::cerr << "[DPU] allocating " << dpu_n << " DPUs" << std::endl;
    auto DPUs = DpuSet::allocate(dpu_n);
    Read_packet rp = get();

    // pipeline management
    while(true) {
        DPUs.copy("my_var", rp);
    }
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
    Read_packet seqs;
    while (true) {
        l = kseq_read(seq);
        if (l <= 0) {
            break;
        }
        auto n = seq->seq.l;
        for (auto i = 0; i < n; i++) {
            char c = seq->seq.s[i];
            c = ::toupper(c);
            if (c=='U') {
                seq->seq.s[i] = 'T';
            } 
        }
        seqs.push_back(seq->seq.s);
        if(seqs.size() == PACKET_SIZE) {
            RNA_read_pool.push_back(seqs);
            seqs.clear();
        }
    }
    if(seqs.size() != 0) {
        RNA_read_pool.push_back(seqs);
    }
    RNA_read_pool.set_end();
    gzclose(fp);
    fp=0;
} 

void Core::allocate_pipeline_worker() {
    for(int i = 0; i < pw_n +1 ; ++i) {
        if(i == 0) {
            pool->enqueue(&Core::read_rna_file, this, readFile);
        }
        else {
            pool->enqueue(&Pipeline_worker::map, pw[i-1]);
        }
    }
} 