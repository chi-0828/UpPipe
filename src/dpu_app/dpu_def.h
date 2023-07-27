#ifndef DPU_DEF_H
#define DPU_DEF_H

#define DPU_PROGRAM "obj/dpu_app/dpu"
#define PACKET_SIZE 12
#define READ_LEN 150
#define PACKET_CAPACITY (READ_LEN*PACKET_SIZE)

#define MAX_table_n 4000000

typedef struct dpu_result{
    int T[18];
    int len;
    int kmer;
}dpu_result;

#endif