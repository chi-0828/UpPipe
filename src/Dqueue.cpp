#include "Dqueue.h"

Read_packet TDeque::pop_front()  {
    std::lock_guard<std::mutex> lock(theMutex);
    Read_packet returnVal = theDeque.front();
    theDeque.pop_front();
    //std::cerr << "[DBG] get, size = " << theDeque.size() << std::endl;
    return returnVal;
}

void TDeque::push_back(Read_packet &val) {
    std::lock_guard<std::mutex> lock(theMutex);
    theDeque.push_back(val);
    //std::cerr << "[DBG] push, size = " << theDeque.size() << std::endl;
}

bool TDeque::empty() {
    std::lock_guard<std::mutex> lock(theMutex);
    //ßstd::cerr << "[DBG] check, size = " << theDeque.size() << std::endl;
    return theDeque.empty();
}

int TDeque::size() {
    std::lock_guard<std::mutex> lock(theMutex);
    return theDeque.size();
}

