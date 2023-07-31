#include "Kmer.hpp"
#include <bitset>
#include <string>
#include <iostream>

using namespace std;


/*// use:  int2bin(a, buffer, buf_size);
// pre:  buf_size >= 8 and buffer has space for buf_size elements
// post: buffer[0,...,7] is the binary representation of a
void int2bin(uint32_t a, char *buffer, int buf_size) {
  //buffer += (buf_size - 1);

  for (int i = 7; i >= 0; i--) {
      *buffer++ = (a & 1) + '0';
      a >>= 1;
  }
}
*/

static const uint64_t twin_table[256] = {
  0xFF, 0xBF, 0x7F, 0x3F, 0xEF, 0xAF, 0x6F, 0x2F,
  0xDF, 0x9F, 0x5F, 0x1F, 0xCF, 0x8F, 0x4F, 0x0F,
  0xFB, 0xBB, 0x7B, 0x3B, 0xEB, 0xAB, 0x6B, 0x2B,
  0xDB, 0x9B, 0x5B, 0x1B, 0xCB, 0x8B, 0x4B, 0x0B,
  0xF7, 0xB7, 0x77, 0x37, 0xE7, 0xA7, 0x67, 0x27,
  0xD7, 0x97, 0x57, 0x17, 0xC7, 0x87, 0x47, 0x07,
  0xF3, 0xB3, 0x73, 0x33, 0xE3, 0xA3, 0x63, 0x23,
  0xD3, 0x93, 0x53, 0x13, 0xC3, 0x83, 0x43, 0x03,
  0xFE, 0xBE, 0x7E, 0x3E, 0xEE, 0xAE, 0x6E, 0x2E,
  0xDE, 0x9E, 0x5E, 0x1E, 0xCE, 0x8E, 0x4E, 0x0E,
  0xFA, 0xBA, 0x7A, 0x3A, 0xEA, 0xAA, 0x6A, 0x2A,
  0xDA, 0x9A, 0x5A, 0x1A, 0xCA, 0x8A, 0x4A, 0x0A,
  0xF6, 0xB6, 0x76, 0x36, 0xE6, 0xA6, 0x66, 0x26,
  0xD6, 0x96, 0x56, 0x16, 0xC6, 0x86, 0x46, 0x06,
  0xF2, 0xB2, 0x72, 0x32, 0xE2, 0xA2, 0x62, 0x22,
  0xD2, 0x92, 0x52, 0x12, 0xC2, 0x82, 0x42, 0x02,
  0xFD, 0xBD, 0x7D, 0x3D, 0xED, 0xAD, 0x6D, 0x2D,
  0xDD, 0x9D, 0x5D, 0x1D, 0xCD, 0x8D, 0x4D, 0x0D,
  0xF9, 0xB9, 0x79, 0x39, 0xE9, 0xA9, 0x69, 0x29,
  0xD9, 0x99, 0x59, 0x19, 0xC9, 0x89, 0x49, 0x09,
  0xF5, 0xB5, 0x75, 0x35, 0xE5, 0xA5, 0x65, 0x25,
  0xD5, 0x95, 0x55, 0x15, 0xC5, 0x85, 0x45, 0x05,
  0xF1, 0xB1, 0x71, 0x31, 0xE1, 0xA1, 0x61, 0x21,
  0xD1, 0x91, 0x51, 0x11, 0xC1, 0x81, 0x41, 0x01,
  0xFC, 0xBC, 0x7C, 0x3C, 0xEC, 0xAC, 0x6C, 0x2C,
  0xDC, 0x9C, 0x5C, 0x1C, 0xCC, 0x8C, 0x4C, 0x0C,
  0xF8, 0xB8, 0x78, 0x38, 0xE8, 0xA8, 0x68, 0x28,
  0xD8, 0x98, 0x58, 0x18, 0xC8, 0x88, 0x48, 0x08,
  0xF4, 0xB4, 0x74, 0x34, 0xE4, 0xA4, 0x64, 0x24,
  0xD4, 0x94, 0x54, 0x14, 0xC4, 0x84, 0x44, 0x04,
  0xF0, 0xB0, 0x70, 0x30, 0xE0, 0xA0, 0x60, 0x20,
  0xD0, 0x90, 0x50, 0x10, 0xC0, 0x80, 0x40, 0x00
};


/*static const uint8_t base_swap[256] = {
	0x00, 0x40, 0x80, 0xc0, 0x10, 0x50, 0x90, 0xd0,
	0x20, 0x60, 0xa0, 0xe0, 0x30, 0x70, 0xb0, 0xf0,
	0x04, 0x44, 0x84, 0xc4, 0x14, 0x54, 0x94, 0xd4,
	0x24, 0x64, 0xa4, 0xe4, 0x34, 0x74, 0xb4, 0xf4,
	0x08, 0x48, 0x88, 0xc8, 0x18, 0x58, 0x98, 0xd8,
	0x28, 0x68, 0xa8, 0xe8, 0x38, 0x78, 0xb8, 0xf8,
	0x0c, 0x4c, 0x8c, 0xcc, 0x1c, 0x5c, 0x9c, 0xdc,
	0x2c, 0x6c, 0xac, 0xec, 0x3c, 0x7c, 0xbc, 0xfc,
	0x01, 0x41, 0x81, 0xc1, 0x11, 0x51, 0x91, 0xd1,
	0x21, 0x61, 0xa1, 0xe1, 0x31, 0x71, 0xb1, 0xf1,
	0x05, 0x45, 0x85, 0xc5, 0x15, 0x55, 0x95, 0xd5,
	0x25, 0x65, 0xa5, 0xe5, 0x35, 0x75, 0xb5, 0xf5,
	0x09, 0x49, 0x89, 0xc9, 0x19, 0x59, 0x99, 0xd9,
	0x29, 0x69, 0xa9, 0xe9, 0x39, 0x79, 0xb9, 0xf9,
	0x0d, 0x4d, 0x8d, 0xcd, 0x1d, 0x5d, 0x9d, 0xdd,
	0x2d, 0x6d, 0xad, 0xed, 0x3d, 0x7d, 0xbd, 0xfd,
	0x02, 0x42, 0x82, 0xc2, 0x12, 0x52, 0x92, 0xd2,
	0x22, 0x62, 0xa2, 0xe2, 0x32, 0x72, 0xb2, 0xf2,
	0x06, 0x46, 0x86, 0xc6, 0x16, 0x56, 0x96, 0xd6,
	0x26, 0x66, 0xa6, 0xe6, 0x36, 0x76, 0xb6, 0xf6,
	0x0a, 0x4a, 0x8a, 0xca, 0x1a, 0x5a, 0x9a, 0xda,
	0x2a, 0x6a, 0xaa, 0xea, 0x3a, 0x7a, 0xba, 0xfa,
	0x0e, 0x4e, 0x8e, 0xce, 0x1e, 0x5e, 0x9e, 0xde,
	0x2e, 0x6e, 0xae, 0xee, 0x3e, 0x7e, 0xbe, 0xfe,
	0x03, 0x43, 0x83, 0xc3, 0x13, 0x53, 0x93, 0xd3,
	0x23, 0x63, 0xa3, 0xe3, 0x33, 0x73, 0xb3, 0xf3,
	0x07, 0x47, 0x87, 0xc7, 0x17, 0x57, 0x97, 0xd7,
	0x27, 0x67, 0xa7, 0xe7, 0x37, 0x77, 0xb7, 0xf7,
	0x0b, 0x4b, 0x8b, 0xcb, 0x1b, 0x5b, 0x9b, 0xdb,
	0x2b, 0x6b, 0xab, 0xeb, 0x3b, 0x7b, 0xbb, 0xfb,
	0x0f, 0x4f, 0x8f, 0xcf, 0x1f, 0x5f, 0x9f, 0xdf,
	0x2f, 0x6f, 0xaf, 0xef, 0x3f, 0x7f, 0xbf, 0xff
};

*/
// use:  km = Kmer();
// pre:
// post: the DNA string in km is AA....AAA (k times A)
Kmer::Kmer() {
    longs = 0UL;
}


// use:  _km = Kmer(km);
// pre:  s[0],...,s[k] are all equal to 'A','C','G' or 'T'
// post: the DNA string in _km and is the same as in km
Kmer::Kmer(const Kmer& o) {
  //memcpy(bytes,o.bytes,MAX_K/4);
  longs = o.longs;
}


// use:  km = Kmer(s);
// pre:  s[0],...,s[k] are all equal to 'A','C','G' or 'T'
// post: the DNA string in km is now the same as s
Kmer::Kmer(const char *s) {
  set_kmer(s);
}


// use:  _km = km;
// pre:
// post: the DNA string in _km and is the same as in km
Kmer& Kmer::operator=(const Kmer& o) {
  if (this != &o) {
    longs = o.longs;
    //memcpy(bytes, o.bytes, MAX_K/4);
  }
  return *this;
}



// use:  km = Kmer();
// pre:
// post: The last 2 bits in the bit array which stores the DNA string have been set to 11
//       which indicates that the km is deleted
void Kmer::set_deleted() {
  longs = DELETE_KMER;
}

// use:  km = Kmer();
// pre:
// post: The last 2 bits in the bit array which stores the DNA string hav been set to 01
//       which indicates that the km is invalid
void Kmer::set_empty() {
  longs = EMPTY_KMER;
}


// use:  b = (km1 < km2);
// pre:
// post: b is true <==> the DNA strings in km1 is alphabetically smaller than
//                      the DNA string in km2
bool Kmer::operator<(const Kmer& o) const {
  if (longs < o.longs) {
    return true;
  }
  if (longs > o.longs) {
    return false;
  }
  return false;
}


// use:  km.set_kmer(s);
// pre:  s[0],...,s[k-1] are all 'A','C','G' or 'T'
// post: The DNA string in km is now equal to s
void Kmer::set_kmer(const char *s)  {
  longs = 0X0UL;
  // std::cerr << "set kmer: ";
  for (size_t i = 0; i < k; ++i) {
    assert(*s != '\0');
    switch(*s) {
      case 'A': longs |= (0x00UL << (2*(k-i-1))); break;
      case 'T': longs |= (0x01UL << (2*(k-i-1))); break;
      case 'C': longs |= (0x02UL << (2*(k-i-1))); break;
      case 'G': longs |= (0x03UL << (2*(k-i-1))); break;
    }
    s++;
  }
  // std::cerr << "\n";
  // std::cerr << longs << "\n";
}

uint64_t Kmer::tobinary() const {
  return longs;
}

// use:  i = km.hash();
// pre:
// post: i is the hash value of km
uint64_t Kmer::hash() const {
  uint64_t ret;
  MurmurHash3_x64_64((const void *)&longs,8,0,&ret);
  return ret;
}


// use:  fw = km.forwardBase(c)
// pre:
// post: fw is the forward kmer from km with last character c,
//       i.e. if the DNA string in km is 'ACGT' and c equals 'T' then
//       the DNA string in fw is 'CGTT'
Kmer Kmer::forwardBase(const char b) const {
  Kmer km(*this);
  uint64_t mask = EMPTY_KMER;
  mask &= ((2UL<<(2*k)) -1);
  
  km.longs = km.longs << 2;
  switch(b) {
      case 'A': km.longs |= (0x00UL); break;
      case 'T': km.longs |= (0x01UL); break;
      case 'C': km.longs |= (0x02UL); break;
      case 'G': km.longs |= (0x03UL); break;
    }
  km.longs &= mask;
  return km;
}


// use:  bw = km.backwardBase(c)
// pre:
// post: bw is the backward kmer from km with first character c,
//       i.e. if the DNA string in km is 'ACGT' and c equals 'T' then
//       the DNA string in bw is 'TACG'
Kmer Kmer::backwardBase(const char b) const {
  Kmer km(*this);
  uint64_t mask = EMPTY_KMER;
  mask &= ((2UL<<(2*k)) -1);
  km.longs = km.longs >> 2;
  switch(b) {
      case 'A': km.longs |= (0x00UL << 2*(k-1)); break;
      case 'T': km.longs |= (0x01UL << 2*(k-1)); break;
      case 'C': km.longs |= (0x02UL << 2*(k-1)); break;
      case 'G': km.longs |= (0x03UL << 2*(k-1)); break;
    }
  km.longs &= mask;
  return km;
}


// use:  km.printBinary();
// pre:
// post: The bits in the binary representation of the
//       DNA string for km has been printed to stdout
std::string Kmer::getBinary() const {

  size_t nlongs = MAX_K/32;
  std::string r;
  r.reserve(64*nlongs);
  for (size_t i = 0; i < nlongs; i++) {
    r.append(std::bitset<64>(longs).to_string<char,std::char_traits<char>,std::allocator<char>>());
  }
  return r;
  /*
  for (size_t i = 0; i < Kmer::k_bytes; i++) {
    int2bin(bytes[i],buff,8);
    printf("%s",buff);
  }

  printf("\n");
  */
}


// use:  km.toString(s);
// pre:  s has space for k+1 elements
// post: s[0,...,k-1] is the DNA string for the Kmer km and s[k] = '\0'
void Kmer::toString(char *s) const {
  size_t j,l;
  for (int i = k - 1; i >= 0; i--) {
    uint64_t kmer = longs;
    kmer = kmer >> (2*i);
    // std::cerr << (kmer& 0x03) << " ";
    switch( kmer & 0x03 ) {
      case 0x00: *s = 'A'; ++s; break;
      case 0x01: *s = 'T'; ++s; break;
      case 0x02: *s = 'C'; ++s; break;
      case 0x03: *s = 'G'; ++s; break;
    }
  }

  *s = '\0';
}

std::string Kmer::toString() const {
  char buf[MAX_K];
  toString(buf);
  return std::string(buf);
}



// use:  set_k(k);
// pre:  this method has not been called before and 0 < k < MAX_K
// post: The Kmer size has been set to k
void Kmer::set_k(unsigned int _k) {
  if(_k == k) {
    return; // ok to call more than once
  }
  assert(_k < MAX_K);
  assert(_k > 0);
  assert(k_bytes == 0); // we can only call this once
  k = _k;
  k_bytes = (_k+3)/4;
  //  k_longs = (_k+15)/16;
  k_modmask = (1 << (2*((k%4)?k%4:4)) )-1;
}


unsigned int Kmer::k = 0;
unsigned int Kmer::k_bytes = 0;
//unsigned int Kmer::k_longs = 0;
unsigned int Kmer::k_modmask = 0;
