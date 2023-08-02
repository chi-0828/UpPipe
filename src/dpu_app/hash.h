#ifndef DPU_HASH_H
#define DPU_HASH_H

#include <stdlib.h>
#include <stdint.h> 
#include "dpu.h"
#include "dpu_def.h"

#define BIG_CONSTANT(x) (x##LLU)

uint64_t getblock ( const uint64_t * p );
//-----------------------------------------------------------------------------
// MurmurHash2, 64-bit versions, by Austin Appleby

// The same caveats as 32-bit MurmurHash2 apply here - beware of alignment 
// and endian-ness issues if used across multiple platforms.

// 64-bit hash for 64-bit platforms

uint64_t murmurhash ( const void * key, int len, uint64_t seed );

uint64_t hash(Kmer* key);


#endif