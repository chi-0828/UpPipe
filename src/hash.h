// source code: https://github.com/explosion/murmurhash/blob/master/murmurhash/MurmurHash2.cpp
//-----------------------------------------------------------------------------
// MurmurHash2 was written by Austin Appleby, and is placed in the public
// domain. The author hereby disclaims copyright to this source code.

// Note - This code makes a few assumptions about how your machine behaves -

// 1. We can read a 4-byte value from any address without crashing
// 2. sizeof(int) == 4

// And it has a few limitations -

// 1. It will not work incrementally.
// 2. It will not produce the same results on little-endian and big-endian
//    machines.
#ifndef HASH_H
#define HASH_H

#include <stdint.h>

//-----------------------------------------------------------------------------
// Block read - on little-endian machines this is a single load,
// while on big-endian or unknown machines the byte accesses should
// still get optimized into the most efficient instruction.
#define BIG_CONSTANT(x) (x##LLU)

static uint64_t getblock ( const uint64_t * p )
{
#if defined(__BYTE_ORDER__) && (__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__)
  return *p;
#else
  const uint8_t *c = (const uint8_t *)p;
  return (uint64_t)c[0] |
	 (uint64_t)c[1] <<  8 |
	 (uint64_t)c[2] << 16 |
	 (uint64_t)c[3] << 24 |
	 (uint64_t)c[4] << 32 |
	 (uint64_t)c[5] << 40 |
	 (uint64_t)c[6] << 48 |
	 (uint64_t)c[7] << 56;
#endif
}
//-----------------------------------------------------------------------------
// MurmurHash2, 64-bit versions, by Austin Appleby

// The same caveats as 32-bit MurmurHash2 apply here - beware of alignment 
// and endian-ness issues if used across multiple platforms.

// 64-bit hash for 64-bit platforms

uint64_t MurmurHash64A ( const void * key, int len, uint64_t seed );

//-----------------------------------------------------------------------------
#endif 