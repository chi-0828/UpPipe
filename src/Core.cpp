#include "Core.h"
#include "Dqueue.h"
#include <chrono>

#ifndef KSEQ_INIT_READY
#define KSEQ_INIT_READY
KSEQ_INIT(gzFile, gzread)
#endif

#include <dpu>
using namespace dpu;

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
std::map<std::set<int16_t>, int> similarity_class;
void insert(Partial_result_host &pr) {
    for(int j = 0; j < PACKET_SIZE; j++) {
        if(pr.at(j).T.size() == 0)
            continue;
        // std::sort(key.begin(), key.end());
        if(similarity_class.find(pr.at(j).T) == similarity_class.end()) {
            similarity_class.insert(std::make_pair(pr.at(j).T, 1));
        }
        else {
            similarity_class.at(pr.at(j).T) ++;
        }
    }
}

// pipeline worker initialization
void Pipeline_worker::operator()(int n, KmerIndex *ki) {
    dpu_n = n;
    KI = ki;
    transfer_to_dpu_buf = new std::deque<Read_packet>(dpu_n, Empty_read_packet);
    dpu_result tmp;
    Partial_result pr_tmp(PACKET_SIZE, tmp);
    transfer_from_dpu_buf = new std::vector<Partial_result>(dpu_n, pr_tmp);
    dpu_result_host tmp_h;
    Partial_result_host pr_tmp_h(PACKET_SIZE, tmp_h);
    partial_result_buf = new std::deque<Partial_result_host>(dpu_n, pr_tmp_h);
    map();
}

// this function keeps the best partial rsults
Partial_result_host Pipeline_worker::compare() {
    for(int i = 0; i < dpu_n; i++) {
        for(int j = 0; j < PACKET_SIZE; j++) {
            if(partial_result_buf->at(i).at(j).kmer < transfer_from_dpu_buf->at(i).at(j).kmer) {
                partial_result_buf->at(i).at(j).kmer = transfer_from_dpu_buf->at(i).at(j).kmer;
                partial_result_buf->at(i).at(j).len = transfer_from_dpu_buf->at(i).at(j).len;
                int16_t *ptr = transfer_from_dpu_buf->at(i).at(j).T;
                partial_result_buf->at(i).at(j).T = std::set<int16_t>(ptr, ptr + transfer_from_dpu_buf->at(i).at(j).len);
            }
            else if(partial_result_buf->at(i).at(j).kmer == transfer_from_dpu_buf->at(i).at(j).kmer) {
                for(int l2 = 0; l2 < transfer_from_dpu_buf->at(i).at(j).len; l2 ++)
                    partial_result_buf->at(i).at(j).T.insert(transfer_from_dpu_buf->at(i).at(j).T[l2]);
                partial_result_buf->at(i).at(j).len = partial_result_buf->at(i).at(j).T.size();
            }
        }
    }
    Partial_result_host final_result = partial_result_buf->back();
    partial_result_buf->pop_back();
    dpu_result_host tmp;
    tmp.kmer = 0;
    Partial_result_host pr_tmp(PACKET_SIZE, tmp);
    partial_result_buf->push_front(pr_tmp);
    return final_result;
}

// start pipeline mapping
void Pipeline_worker::map() {
    auto DPUs = DpuSet::allocate(dpu_n);
    auto dpu = DPUs.dpus();
    DPUs.load(DPU_PROGRAM);

    // transfer database and parameters
    DPUs.copy("first_tid" ,KI->first_tid_buf);
    DPUs.copy("t_max" ,KI->t_max_buf);
    DPUs.copy("k" ,KI->k_buf);
    DPUs.copy("table" ,KI->table_buf);
    DPUs.copy("size_" ,KI->size_buf);

    // pipeline management
    int stage = 0;
    bool end = false;
    double get_read_time = 0, CPU_DPU_time = 0, DPU_run_time = 0, DPU_CPU_time = 0, compare_time = 0, insert_time = 0;
    while(stage < dpu_n) {
        // ====== start CPU-DPU transfer =====
        auto t_start = std::chrono::high_resolution_clock::now();
        Read_packet rp = get(end);
        if(end)
            stage ++;
        transfer_to_dpu_buf->pop_back();
        transfer_to_dpu_buf->push_front(rp);
        std::vector<Read_packet> to_dpu_buf({transfer_to_dpu_buf->begin(), transfer_to_dpu_buf->end()});
        auto t_end = std::chrono::high_resolution_clock::now();
        get_read_time += std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
        
        t_start = std::chrono::high_resolution_clock::now();
        DPUs.copy("reads", to_dpu_buf);
        t_end = std::chrono::high_resolution_clock::now();
        CPU_DPU_time += std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
        // ====== end CPU-DPU transfer =====
        
        // ====== start DPU execution =====
        t_start = std::chrono::high_resolution_clock::now();
        DPUs.exec();
        t_end = std::chrono::high_resolution_clock::now();
        DPU_run_time += std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
        // ====== end DPU execution =====

        // to print DPU-printf
        // DPUs.log(std::cerr);
        
        // ====== start DPU-CPU transfer =====
        t_start = std::chrono::high_resolution_clock::now();
        DPUs.copy(*transfer_from_dpu_buf, "results");
        t_end = std::chrono::high_resolution_clock::now();
        DPU_CPU_time += std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
        

        t_start = std::chrono::high_resolution_clock::now();
        Partial_result_host final_result = compare();
        t_end = std::chrono::high_resolution_clock::now();
        compare_time += std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
        

        t_start = std::chrono::high_resolution_clock::now();
        insert(final_result);
        t_end = std::chrono::high_resolution_clock::now();
        insert_time += std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
        // ====== end DPU-CPU transfer =====
    }
    // uncomment here to output breakdown execution time 
    // std::cerr << "get reads time " << get_read_time << " s\n";
    // std::cerr << "CPU-DPU time " << CPU_DPU_time << " s\n";
    // std::cerr << "DPU run time " << DPU_run_time << "s\n";
    // std::cerr << "DPU-CPU time " << DPU_CPU_time << "s\n";
    // std::cerr << "compare result time " << compare_time << "s\n";
    // std::cerr << "insert result time " << insert_time << "s\n";
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
    if(!fp) {
        std::cerr << "open file error\n"; 
        exit(1);
    }
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
    std::cerr << "[UpPipe] start processing ... \n";
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
    std::cerr << "[output] write similarity class ... " ;
    std::ofstream out;
    out.open(outputFile, std::ios::out);
    if (!out.is_open()) {
		std::cerr << "Error: output file \"" << outputFile << "\" could not be opened!\n";
		exit(1);
	}
    for(auto& entry: similarity_class) {
        // class
        for(auto& t: entry.first)
            out << t << " ";
        out << ": ";
        // count
        out << entry.second << "\n";
    }
    std::cerr << "done\n";
} 