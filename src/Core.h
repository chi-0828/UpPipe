#ifndef CORE_H
#define CORE_H

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <thread>
#include <mutex>
#include <algorithm>
#include <ctype.h>
#include <zlib.h>
#include "KmerIndex.h"
#include <omp.h>
#include "kseq.h"
#include "Dqueue.h"
#include "dpu_app/dpu_def.h"

#define TIME_PERF 0
#define Empty_read_packet Read_packet(PACKET_CAPACITY, 'N')

// global similarity class record 
class Global_Record {
public:
    std::map<std::set<int16_t>, uint32_t> similarity_class;
	std::mutex mtx;  

	void insert(Partial_result_host &pr);
};

class Pipeline_worker {
public:
    Pipeline_worker(int n, KmerIndex *ki, Global_Record *gr, TDeque *rb, bool *end, int pid) : 
        dpu_n(n), KI(ki), record(gr), read_buffer(rb), input_end(end), pid(pid) {
            transfer_to_dpu_buf = new std::deque<Read_packet>(dpu_n, Empty_read_packet);
            dpu_result tmp;
            Partial_result pr_tmp(PACKET_SIZE, tmp);
            transfer_from_dpu_buf = new std::vector<Partial_result>(dpu_n, pr_tmp);
            dpu_result_host tmp_h;
            Partial_result_host pr_tmp_h(PACKET_SIZE, tmp_h);
            partial_result_buf = new std::deque<Partial_result_host>(dpu_n, pr_tmp_h);
            map();
        }

    ~Pipeline_worker() {
        delete transfer_to_dpu_buf;
        delete transfer_from_dpu_buf;
        delete partial_result_buf;
    }
    int pid;
    KmerIndex *KI;
	Global_Record *record;
    TDeque *read_buffer;
    // std::atomic<bool> *input_end;
    bool *input_end;
    int dpu_n;
    void operator()(int n, KmerIndex *ki, Global_Record *gr, TDeque *rb, bool *end);
    void map();
    Read_packet get_read_packet(bool &sys_end);
    Partial_result_host compare();

    std::deque<Read_packet> *transfer_to_dpu_buf;
    std::vector<Partial_result> *transfer_from_dpu_buf;
    std::deque<Partial_result_host> *partial_result_buf;
};

class Core {
public:
    Core(int pw_n, int dpu_n, KmerIndex *KI, std::string &readFile, std::string &outputFile) : pw_n(pw_n), dpu_n(dpu_n), KI(KI), readFile(readFile), outputFile(outputFile) {
        std::cerr << "[UpPipe] start alignment by k = " << KI->k << std::endl;
		gr = new Global_Record();
        read_buffers = new std::vector<TDeque>(pw_n);
        allocate_pipeline_worker();
    }
    ~Core() {
		delete gr;
        delete read_buffers;
        std::cerr << "[UpPipe] end alignment" << std::endl;
    }

    int pw_n; // number of pipeline workers
    int dpu_n; // number of DPUs in a pipeline worker
    KmerIndex *KI;
    std::string readFile;
    std::string outputFile;
	Global_Record *gr;

    std::vector<TDeque> *read_buffers;
    // std::atomic<bool> input_end;
    bool input_end;

    void read_rna_file();
    void write_out_file();
    void allocate_pipeline_worker();
};

#endif 