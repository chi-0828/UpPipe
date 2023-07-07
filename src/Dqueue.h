#ifndef DQUEUE_H
#define DQUEUE_H

#include "Core.h"

class TDeque {
public:
    Read_packet pop_front();

    void push_back(Read_packet &val);

    bool empty();

    void set_end();
    bool get_end();

    std::deque<Read_packet> theDeque;
    std::mutex theMutex;
    std::mutex endMutex;
    bool end = false;
};

#endif 