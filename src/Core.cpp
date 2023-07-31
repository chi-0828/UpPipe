#include "Core.h"
#include "Dqueue.h"
#include <chrono>

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

// pipeline worker 
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

Partial_result_host Pipeline_worker::compare() {
    for(int i = 0; i < dpu_n; i++) {
        for(int j = 0; j < PACKET_SIZE; j++) {
            assert(transfer_from_dpu_buf->at(i).at(j).kmer >= 0);
            if(transfer_from_dpu_buf->at(i).at(j).len < 0 || transfer_from_dpu_buf->at(i).at(j).len > T_LEN) {
                std::cerr << transfer_from_dpu_buf->at(i).at(j).len << " < 0 !\n";
                for(int ll = 0; ll < T_LEN; ll ++)
                    std::cerr << transfer_from_dpu_buf->at(i).at(j).T[ll] << " ";
                std::cerr << "\n";
            }
            if(partial_result_buf->at(i).at(j).kmer < transfer_from_dpu_buf->at(i).at(j).kmer) {
                partial_result_buf->at(i).at(j).kmer = transfer_from_dpu_buf->at(i).at(j).kmer;
                partial_result_buf->at(i).at(j).len = transfer_from_dpu_buf->at(i).at(j).len;
                int16_t *ptr = transfer_from_dpu_buf->at(i).at(j).T;
                partial_result_buf->at(i).at(j).T = std::set<int16_t>(ptr, ptr + transfer_from_dpu_buf->at(i).at(j).len);
            }
            else if(partial_result_buf->at(i).at(j).kmer == transfer_from_dpu_buf->at(i).at(j).kmer) {
                for(int l2 = 0; l2 < transfer_from_dpu_buf->at(i).at(j).len; l2 ++)
                    partial_result_buf->at(i).at(j).T.insert(transfer_from_dpu_buf->at(i).at(j).T[l2]);
                // for(auto& map_item: check_map)
                //     pushed_items.push_back(map_item);
                partial_result_buf->at(i).at(j).len = partial_result_buf->at(i).at(j).T.size();
            }
            //assert(partial_result_buf->at(i).at(j).len >= 0);
            //assert(partial_result_buf->at(i).at(j).len <= T_LEN);
        }
    }
    Partial_result_host final_result = partial_result_buf->front();
    partial_result_buf->pop_front();
    dpu_result_host tmp;
    memset(&tmp, 0, sizeof(dpu_result_host));
    Partial_result_host pr_tmp(PACKET_SIZE, tmp);
    partial_result_buf->push_back(pr_tmp);
    return final_result;
}

void Pipeline_worker::map() {
    // std::cerr << "[DPU] allocating " << dpu_n << " DPUs" << std::endl;
    auto DPUs = DpuSet::allocate(dpu_n);
    auto dpu = DPUs.dpus();
    // std::cerr << "[DPU] loading " << DPU_PROGRAM << std::endl;
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
    
    double get_read_time = 0, CPU_DPU_time = 0, DPU_run_time = 0, DPU_CPU_time = 0, compare_time = 0, insert_time = 0;
    while(stage < dpu_n) {
        auto t_start = std::chrono::high_resolution_clock::now();
        // std::cerr << "[pipe] stage " << stage << " end " << end << std::endl;
        Read_packet rp = get(end);
        if(end)
            stage ++;
        transfer_to_dpu_buf->pop_front();
        transfer_to_dpu_buf->push_back(rp);
        std::vector<Read_packet> to_dpu_buf({transfer_to_dpu_buf->begin(), transfer_to_dpu_buf->end()});
        auto t_end = std::chrono::high_resolution_clock::now();
        get_read_time += std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
        
        // std::cerr << "[pipe] copy to " << std::endl;
        t_start = std::chrono::high_resolution_clock::now();
        DPUs.copy("reads", to_dpu_buf);
        t_end = std::chrono::high_resolution_clock::now();
        CPU_DPU_time += std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
        
        // std::cerr << "+===============\n";
        t_start = std::chrono::high_resolution_clock::now();
        DPUs.exec();
        t_end = std::chrono::high_resolution_clock::now();
        DPU_run_time += std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
       
        // DPUs.log(std::cerr);
        t_start = std::chrono::high_resolution_clock::now();
        DPUs.copy(*transfer_from_dpu_buf, "result");
        t_end = std::chrono::high_resolution_clock::now();
        DPU_CPU_time += std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
        
        // std::cerr << "-===============\n";

        t_start = std::chrono::high_resolution_clock::now();
        Partial_result_host final_result = compare();
        t_end = std::chrono::high_resolution_clock::now();
        compare_time += std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
        

        t_start = std::chrono::high_resolution_clock::now();
        insert(final_result);
        t_end = std::chrono::high_resolution_clock::now();
        insert_time += std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
    }
    // std::cerr << "get reads time " << get_read_time << " s\n";
    // std::cerr << "CPU-DPU time " << CPU_DPU_time << " s\n";
    // std::cerr << "DPU run time " << DPU_run_time << "\n";
    // std::cerr << "DPU-CPU time " << DPU_CPU_time << "\n";
    // std::cerr << "compare result time " << compare_time << "\n";
    // std::cerr << "insert result time " << insert_time << "\n";
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
    std::cerr << "[output] write similarity class ...\n" ;
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
        out << " = ";
        // count
        out << entry.second << "\n";
    }
    out << "\n";
} 