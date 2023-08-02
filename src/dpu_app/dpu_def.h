#ifndef DPU_DEF_H
#define DPU_DEF_H
#include <stdint.h>
#define DPU_PROGRAM "obj/dpu_app/dpu"
#define PACKET_SIZE 12
#define READ_LEN 150
#define PACKET_CAPACITY (READ_LEN*PACKET_SIZE)
#define T_LEN 104
#define MAX_table_n 5000000

typedef struct dpu_result{
    int32_t kmer;
    int32_t len;
    int16_t T[T_LEN];
}dpu_result;

#endif