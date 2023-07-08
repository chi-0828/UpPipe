#ifndef DPU_DEF_H
#define DPU_DEF_H

#define DPU_PROGRAM "obj/dpu_app/dpu"
#define PACKET_SIZE 12
#define MAX_READ_LEN 150
#define MAX_READ_LEN 150
#define MAX_PACKET_SIZE (MAX_READ_LEN*PACKET_SIZE)

#define MAX_table_n 4000000

typedef struct dpu_args {
    // k-mer size
    int k;
    // mapping or intersection
    int run;
}dpu_args;

typedef struct dpu_result{
    int T[18];
    int len;
    int kmer;
}dpu_result;

#endif