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
bool end = false;
TDeque RNA_read_pool;
Read_packet get(bool &end) {
    if(end)
        return Read_packet(0);

    Read_packet rp;

    while(RNA_read_pool.empty() && !RNA_read_pool.get_end()) ;
    
    s
    if(RNA_read_pool.get_end()) {
        rp = Read_packet(0);
        end = true;
    }
    else
        rp = RNA_read_pool.pop_front();

    return rp;
}

// global similarity class record
std::map<std::vector<int>, int> similarity_class;
void insert(Partial_result &pr) {
    for(int j = 0; j < PACKET_SIZE; j++) {
        std::vector<int> key(pr.at(j).T, pr.at(j).T + pr.at(j).len);
        if(similarity_class.find(key) == similarity_class.end()) {
            similarity_class.insert(std::make_pair(key, 1));
        }
        else {
            similarity_class.at(key) ++;
        }
    }
}

// pipeline worker 
Partial_result Pipeline_worker::compare() {
    for(int i = 0; i < dpu_n; i++) {
        for(int j = 0; j < PACKET_SIZE; j++) {
            if(partial_result_buf->at(i).at(j).kmer < transfer_from_dpu_buf->at(i).at(j).kmer) {
                partial_result_buf->at(i).at(j) = transfer_from_dpu_buf->at(i).at(j);
            }
        }
    }
    Partial_result final_result = partial_result_buf->front();
    partial_result_buf->pop_front();
    dpu_result tmp;
    tmp.kmer = 0;
    Partial_result pr_tmp(PACKET_SIZE, tmp);
    partial_result_buf->push_back(pr_tmp);
    return final_result;
}

void Pipeline_worker::map() {
    std::cerr << "[DPU] allocating " << dpu_n << " DPUs" << std::endl;
    auto DPUs = DpuSet::allocate(dpu_n);

    // pipeline management
    int stage = 0;
    bool end = false;
    while(stage <= dpu_n) {
        if(end)
            stage ++;
        Read_packet rp = get(end);
        transfer_to_dpu_buf->pop_front();
        transfer_to_dpu_buf->push_back(rp);
        std::vector<Read_packet> to_dpu_buf = std::vector<Read_packet>({transfer_to_dpu_buf->begin(), transfer_to_dpu_buf->end()});
        DPUs.copy("Read_packet", to_dpu_buf);
        DPUs.exec();
        DPUs.copy(*transfer_from_dpu_buf, "Partial_result");

        Partial_result final_result = compare();
        insert(final_result);
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
    // output similarity class
} 