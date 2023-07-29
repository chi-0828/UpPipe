#ifndef DPU_DEF_H
#define DPU_DEF_H
#include <stdint.h>
#define DPU_PROGRAM "obj/dpu_app/dpu"
#define PACKET_SIZE 12
#define READ_LEN 8
#define PACKET_CAPACITY (READ_LEN*PACKET_SIZE)

#define MAX_table_n 4000000

typedef struct dpu_result{
    int16_t T[31];
    int16_t len;
}dpu_result;

#endif