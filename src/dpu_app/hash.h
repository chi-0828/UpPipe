#ifndef HASH_H
#define HASH_H

#include <stdlib.h>
#include <stdint.h> 
#include "dpu.h"
#include "dpu_def.h"

#undef get16bits
#if (defined(__GNUC__) && defined(__i386__)) || defined(__WATCOMC__) \
  || defined(_MSC_VER) || defined (__BORLANDC__) || defined (__TURBOC__)
#define get16bits(d) (*((const uint16_t *) (d)))
#endif

#if !defined (get16bits)
#define get16bits(d) ((((uint32_t)(((const uint8_t *)(d))[1])) << 8)\
                       +(uint32_t)(((const uint8_t *)(d))[0]) )
#endif

uint32_t SuperFastHash (const char *data, int len);

//void MurmurHash3_x64_32 ( const void * key, int len, uint32_t seed, void * out );
void MurmurHash3_x64_64 ( const void *key, int len, uint32_t seed, void *out );

// hash table
typedef struct KmerHashTable {
    // int *table_int;
    // Kmer *table_kmer;
    Kmer empty;
    size_tt size_;
    __mram_ptr uint64_t *table_ptr;
}KmerHashTable;

uint64_t hash(Kmer* key);

#endif