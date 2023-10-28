#ifndef DQUEUE_H
#define DQUEUE_H

#include <atomic>
#include <vector>
#include <set>
#include <deque>
#include <mutex>
#include "dpu_app/dpu_def.h"

// re-define dpu_result in host for larger memory size
// original dpu_result is defined in dpu_app/dpu_def.h
typedef struct dpu_result_host{
    int32_t kmer;
    int32_t len;
    std::set<int16_t> T;
}dpu_result_host;

using Read_packet = std::vector<char>;
using Partial_result =  std::vector<dpu_result>;
using Partial_result_host =  std::vector<dpu_result_host>;

class TDeque {
public:
    Read_packet pop_front();

    void push_back(Read_packet &val);

    bool empty();

    int size();

    std::deque<Read_packet> theDeque;
    std::mutex theMutex;
};


#endif 