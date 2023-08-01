#include "dpu.h"
#include <mram.h>
#include <assert.h>
#include "mutex.h"
#include <defs.h>
#include "dpu_def.h"

extern int32_t k;

void KmerIteratornew(KmerIterator* iterator, const char *s){
    iterator->s_ = s;
    iterator->p_ = 0;
    iterator->invalid_ = 0;
    find_next(iterator,-1,-1,0);
}

void KmerIteratornew_NULL(KmerIterator* iterator){
    iterator->s_ = NULL;
    iterator->invalid_ = 1;
}


void kmernew(Kmer* kmer, const char *s){
    set_kmer(s, kmer);
}
void kmernew_copy(Kmer* kmer, Kmer* o){
    for (size_tt i = 0; i < MAX_K/32; i++) {
        kmer->longs[i] = o->longs[i];
    }
}
void Kmerfree(Kmer* kmer){
  buddy_free(kmer);
}

void KmerIterator_move(KmerIterator* kmerIterator) {
    int pos_ = kmerIterator->p_;
    if (!kmerIterator->invalid_) {
        if (kmerIterator->s_[pos_+k] == 0) {
            kmerIterator->invalid_ = 1;
        } else {
            find_next(kmerIterator, pos_,pos_+k-1,1);
        }
    }
}

void find_next(KmerIterator* kmerIterator, size_tt i, size_tt j, uint8_t last_valid) {
  ++i;
  ++j;

  while (kmerIterator->s_[j] != 0) {
    char c = kmerIterator->s_[j] & 0xDF; // mask lowercase bit
    if (c == 'A' || c == 'C' || c == 'G' || c == 'T') {
      if (last_valid) {
        forwardBase(&(kmerIterator->kmer), c);
        break; // default case,
      } else {
        if (i + k - 1 == j) {
          Kmer k_tmp;
          kmernew(&k_tmp, kmerIterator->s_+i);
          kmernew_copy(&kmerIterator->kmer, &k_tmp);
          last_valid = 1;
          break; // create k-mer from scratch
        } else {
          ++j;
        }
      }
    } else {
      ++j;
      i = j;
      last_valid = 0;
    }
  }
  if (i+k-1 == j && kmerIterator->s_[j] != 0) {
    kmerIterator->p_ = i;
  } else {
    kmerIterator->invalid_ = 1;
  }
  //toString(kmerIterator->kmer);
}

void forwardBase(Kmer* kmer, const char b) {
  uint64_t mask = EMPTY_KMER;
  mask &= ((1UL<<(2*k)) -1);
	kmer->longs[0] = kmer->longs[0] << 2;
	switch(b) {
		case 'A': kmer->longs[0] |= (0x00UL); break;
		case 'T': kmer->longs[0] |= (0x01UL); break;
		case 'C': kmer->longs[0] |= (0x02UL); break;
		case 'G': kmer->longs[0] |= (0x03UL); break;
	}
  kmer->longs[0] &= mask;
	// std::cerr << kmer->longs[0]getBinary() << "\n" ;
}
// < 
uint8_t Kmer_cmp(const Kmer* a, const Kmer* b){
  for (size_tt i = 0; i < MAX_K/32; ++i) {
    if (a->longs[i] < b->longs[i]) {
      return 1;
    }
    if (a->longs[i] > b->longs[i]) {
      return 0;
    }
  }
  return 0;
}

inline uint8_t Kmer_eq(const Kmer* a, const Kmer* b) {
  for (size_tt i = 0; i < MAX_K/32; i++) {
    if (a->longs[i] != b->longs[i]) {
      return 0;
    }
  }
  return 1;
}

Kmer rep(Kmer km) {
  Kmer tw;
  kmernew_copy(&tw, &km);

  twin(&tw);
  if(Kmer_cmp(&tw, &km)){
    kmernew_copy(&km, &tw);
  }
  return km;
}

void twin(Kmer* km) {
  size_tt nlongs = (k+31)/32;
  // may have bug
  uint64_t this_longs[1];
  for (size_tt i = 0; i < nlongs; i++){
    this_longs[i] = km->longs[i];
  }
  for (size_tt i = 0; i < nlongs; i++) {
    uint64_t v = this_longs[i];
    km->longs[nlongs-1-i] =
      (twin_table[v & 0xFF] << 56) |
      (twin_table[(v>>8) & 0xFF] << 48) |
      (twin_table[(v>>16) & 0xFF] << 40) |
      (twin_table[(v>>24) & 0xFF] << 32) |
      (twin_table[(v>>32) & 0xFF] << 24) |
      (twin_table[(v>>40) & 0xFF] << 16) |
      (twin_table[(v>>48) & 0xFF] << 8)  |
      (twin_table[(v>>56)]);
  }
  size_tt shift = (k%32) ? 2*(32-(k%32)) : 0 ;
  uint64_t shiftmask = (k%32) ? (((1ULL<< shift)-1) << (64-shift)) : 0ULL;
  km->longs[0] = km->longs[0] << shift;

  for (size_tt i = 1; i < nlongs; i++) {
    km->longs[i-1] |= (km->longs[i] & shiftmask) >> (64-shift);
    km->longs[i] = km->longs[i] << shift;
  }
}

void set_kmer(const char *s, Kmer* kmer){
  memset(kmer->bytes,0,MAX_K/4);
  for (size_t i = 0; i < k; ++i) {
    assert(*s != '\0');
    switch(*s) {
      case 'A': kmer->longs[0] |= (0x00UL << (2*(k-i-1))); break;
      case 'T': kmer->longs[0] |= (0x01UL << (2*(k-i-1))); break;
      case 'C': kmer->longs[0] |= (0x02UL << (2*(k-i-1))); break;
      case 'G': kmer->longs[0] |= (0x03UL << (2*(k-i-1))); break;
    }
    s++;
  }
}

uint8_t KmerIterator_cmp(KmerIterator* a, KmerIterator* b){
  if (a->invalid_  || b->invalid_) {
    return a->invalid_ && b->invalid_;
  } else {
    return (a->p_ == b->p_) && (a->s_ == b->s_);
  }
}

void toString(Kmer* kmer) {

  char buf[MAX_K];
  char*s = buf;
  for (int i = k - 1; i >= 0; i--) {
    uint64_t v = kmer->longs[0];
    v = v >> (2*i);
    // std::cerr << (v& 0x03) << " ";
    switch( v & 0x03 ) {
      case 0x00: *s = 'A'; ++s; break;
      case 0x01: *s = 'T'; ++s; break;
      case 0x02: *s = 'C'; ++s; break;
      case 0x03: *s = 'G'; ++s; break;
    }
  }

  *s = '\0';
  printf("%s \n", buf);
}

void set_empty(Kmer* kmer) {
  kmer->longs[0] = EMPTY_KMER;
}



