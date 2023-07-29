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
Read_packet Empty_read_packet(PACKET_CAPACITY, 'N');
Read_packet get(bool &end) {
    if(end)
        return Empty_read_packet;

    Read_packet rp;

    while(RNA_read_pool.empty() && !RNA_read_pool.get_end()) ;
    
    if(RNA_read_pool.get_end() && RNA_read_pool.empty()) {
        rp = Empty_read_packet;
        end = true;
    }
    else if (RNA_read_pool.get_end() && RNA_read_pool.size() <= 1) {
        rp = RNA_read_pool.pop_front();
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
void Pipeline_worker::operator()(int n, KmerIndex *ki) {
    dpu_n = n;
    KI = ki;
    transfer_to_dpu_buf = new std::deque<Read_packet>(dpu_n, Empty_read_packet);
    dpu_result tmp;
    Partial_result pr_tmp(PACKET_SIZE, tmp);
    transfer_from_dpu_buf = new std::vector<Partial_result>(dpu_n, pr_tmp);
    partial_result_buf = new std::deque<Partial_result>(dpu_n, pr_tmp);
    map();
}

Partial_result Pipeline_worker::compare() {
    for(int i = 0; i < dpu_n; i++) {
        for(int j = 0; j < PACKET_SIZE; j++) {
            if(partial_result_buf->at(i).at(j).len < transfer_from_dpu_buf->at(i).at(j).len) {
                partial_result_buf->at(i).at(j) = transfer_from_dpu_buf->at(i).at(j);
            }
        }
    }
    Partial_result final_result = partial_result_buf->front();
    partial_result_buf->pop_front();
    dpu_result tmp;
    memset(&tmp, 0, sizeof(dpu_result));
    Partial_result pr_tmp(PACKET_SIZE, tmp);
    partial_result_buf->push_back(pr_tmp);
    return final_result;
}

void Pipeline_worker::map() {
    std::cerr << "[DPU] allocating " << dpu_n << " DPUs" << std::endl;
    auto DPUs = DpuSet::allocate(dpu_n);
    auto dpu = DPUs.dpus();
    std::cerr << "[DPU] loading " << DPU_PROGRAM << std::endl;
    DPUs.load(DPU_PROGRAM);

    // transfer database and parameters
    DPUs.copy("kmer_max" ,KI->kmer_max_buf);
    DPUs.copy("t_max" ,KI->t_max_buf);
    DPUs.copy("k" ,KI->k_buf);
    DPUs.copy("table" ,KI->table_buf);
    DPUs.copy("size_" ,KI->size_buf);

    // pipeline management
    int stage = 0;
    bool end = false;
    while(stage < dpu_n) {
        // std::cerr << "[pipe] stage " << stage << " end " << end << std::endl;
        Read_packet rp = get(end);
        if(end)
            stage ++;
        transfer_to_dpu_buf->pop_front();
        transfer_to_dpu_buf->push_back(rp);
        std::vector<Read_packet> to_dpu_buf = std::vector<Read_packet>({transfer_to_dpu_buf->begin(), transfer_to_dpu_buf->end()});
        // std::cerr << "[pipe] copy to " << std::endl;
        DPUs.copy("reads", to_dpu_buf);
        std::cerr << "================\n";
        DPUs.exec();
        DPUs.log(std::cerr);
        std::cerr << "================\n";
        DPUs.copy(*transfer_from_dpu_buf, "result");
        // exit(1);

        // Partial_result final_result = compare();
        // insert(final_result);
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
    int rp_size = 0;
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
        seqs.insert(seqs.end(), seq->seq.s, seq->seq.s + READ_LEN);
        rp_size ++;
        if(rp_size == PACKET_SIZE) {
            RNA_read_pool.push_back(seqs);
            seqs.clear();
            rp_size = 0;
        }
    }
    if(rp_size != 0) {
        Read_packet rmp(READ_LEN*(PACKET_SIZE-rp_size), 'N');
        seqs.insert( seqs.end() , rmp.begin() , rmp.end());
        RNA_read_pool.push_back(seqs);
    }
    RNA_read_pool.set_end();
    gzclose(fp);
    fp=0;
} 

void Core::allocate_pipeline_worker() {
    std::vector<std::thread> UpPipe_sys;
    for(int i = 0; i < pw_n +1 ; ++i) {
        if(i == 0) {
            UpPipe_sys.emplace_back(std::thread(&Core::read_rna_file, this, std::ref(readFile)));
        }
        else {
            UpPipe_sys.emplace_back(std::thread(Pipeline_worker(), dpu_n, KI));
        }
    }
    for (std::thread & th : UpPipe_sys) {
        if (th.joinable())
            th.join();
    }
    // output similarity class
} 