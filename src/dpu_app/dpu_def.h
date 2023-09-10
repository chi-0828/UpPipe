#ifndef DPU_DEF_H
#define DPU_DEF_H
#include <stdint.h>
#define DPU_PROGRAM "obj/dpu_app/dpu"
#define PACKET_SIZE 24
#define READ_LEN 101
#define WRAM_READ_LEN 112
#define PACKET_CAPACITY (READ_LEN*PACKET_SIZE)
#define T_LEN 16
#define COUNT_LEN 200
#define MAX_TABLE_N 3000000
#define WRAM_PREFETCH_SIZE 16

typedef struct dpu_result{
    int32_t kmer;
    int32_t len;
    int16_t T[T_LEN];
}dpu_result;

#endif