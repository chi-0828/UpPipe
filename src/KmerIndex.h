#ifndef KALLISTO_KMERINDEX_H
#define KALLISTO_KMERINDEX_H

#ifdef DPU_HEAD
#define DPU_HEAD 1
#else
#include "dpu_app/dpu_def.h"
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <unordered_map>
#include <stdint.h>
#include <ostream>
#include <map>
#include <deque>
#include "common.h"
#include "Kmer.hpp"
#include "KmerIterator.hpp"
#include "KmerHashTable.h"
#include "hash.hpp"

inline int W(int k) {
  if(k < 10)
    return k;
  if(k < 20)
    return k -6;
  return k -16;
}
inline Kmer get_mini(std::deque<Kmer> &dq) {
  Kmer seed;
  uint64_t seed_binary = EMPTY_KMER;
  for(auto& kmer: dq) {
    if(seed_binary >= kmer.tobinary()) {
      seed_binary = kmer.tobinary();
      seed = kmer;
    }
  }
  return seed;
}

std::string revcomp(const std::string s);


struct TRInfo {
  int trid;
  int start;
  int stop; //exclusive [start,stop)
  bool sense; // true for sense, false for anti-sense
};

using EcMap = std::vector<std::vector<int>>; //std::unordered_map<int, std::vector<int>>;

struct SortedVectorHasher {
  size_t operator()(const std::vector<int>& v) const {
    uint64_t r = 0;
    int i=0;
    for (auto x : v) {
      uint64_t t;
      MurmurHash3_x64_64(&x,sizeof(x), 0,&t);
      t = (x>>i) | (x<<(64-i));
      r = r ^ t;
      i = (i+1)%64;
    }
    return r;
  }
};

struct KmerEntry {
  int32_t contig; // id of contig
  uint32_t _pos; // 0-based forward distance to EC-junction
  int32_t contig_length;

  KmerEntry() : contig(-1), _pos(0xFFFFFFF), contig_length(0) {}
  KmerEntry(int id, int length, int pos, bool isFw) : contig(id), contig_length(length) {
    setPos(pos);
    setDir(isFw);
  }

  inline int getPos() const {return (_pos & 0x0FFFFFFF);}
  inline int isFw() const  {return (_pos & 0xF0000000) == 0; }
  inline void setPos(int p) {_pos = (_pos & 0xF0000000) | (p & 0x0FFFFFFF);}
  inline void setDir(bool _isFw) {_pos = (_pos & 0x0FFFFFFF) | ((_isFw) ? 0 : 0xF0000000);}
  inline int getDist(bool fw) const {
    if (isFw() == fw) {
      return (contig_length - 1 - getPos());
    } else {
      return getPos();
    }
  }
};

struct ContigToTranscript {
  int trid;
  int pos; 
  bool sense; // true for sense, 
};

struct Contig {
  int id; // internal id
  int length; // number of k-mers
  int ec;
  std::string seq; // sequence
  std::vector<ContigToTranscript> transcripts;
};

struct DBGraph {
  std::vector<int> ecs; // contig id -> ec-id
  std::vector<Contig> contigs; // contig id -> contig
//  std::vector<pair<int, bool>> edges; // contig id -> edges
};

struct KmerIndex {
  
  KmerIndex(const ProgramOptions& opt) : dpu_n(opt.dpu_n), k(opt.k), num_trans(0), target_seqs_loaded(false), t_max(1), kmer_max(0) {
	  // hash_tables = new std::vector<std::map<Kmer, std::vector<int16_t>>>(dpu_n);
  }
  ~KmerIndex() {
	  // delete hash_tables;
  }

  void hashtable_aligner();
  void Build(const ProgramOptions& opt);
  void Buildhashtable(const std::vector<std::string>& seqs);
  void BuildEquivalenceClasses(const ProgramOptions& opt, const std::vector<std::string>& seqs);


  // output methods
  void write(const std::string& index_out, bool writeKmerTable = true);
  
  // note opt is not const
  // load methods
  void load(ProgramOptions& opt);


  int dpu_n;
  int k; // k-mer size used
  int num_trans; // number of targets

  KmerHashTable<std::set<int16_t>, KmerHash> kmap;
  EcMap ecmap;
  DBGraph dbGraph;
  std::unordered_map<std::vector<int>, int, SortedVectorHasher> ecmapinv;
  const size_t INDEX_VERSION = 10; // increase this every time you change the fileformat

  std::vector<int> target_lens_;

  std::vector<std::string> target_names_;
  std::vector<std::string> target_seqs_; // populated on demand
  bool target_seqs_loaded;

	std::vector<int32_t> kmer_max_buf;
	std::vector<std::vector<int32_t>> t_max_buf;
	std::vector<int32_t> k_buf;
  std::vector<std::vector<size_t>> size_buf;
	std::vector<std::vector<uint64_t>> table_buf;
	// std::vector<std::map<Kmer, std::vector<int16_t>>>* hash_tables;
	int t_max;
	int kmer_max;
};








#endif // KALLISTO_KMERINDEX_H
