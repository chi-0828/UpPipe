#include "KmerIndex.h"
#include <algorithm>
#include <random>
#include <sstream>
#include <ctype.h>
#include <zlib.h>
#include <unordered_set>
#include "kseq.h"

#ifndef KSEQ_INIT_READY
#define KSEQ_INIT_READY
KSEQ_INIT(gzFile, gzread)
#endif

void roundup_4(int32_t& n) {
	n = (n + 3) & ~0x03;
}

// helper functions
// pre: u is sorted
bool isUnique(const std::vector<int>& u) {
  for (int j = 1; j < u.size(); j++) {
    if (u[j-1] == u[j]) {
      return false;
    }
  }
  return true;
}

std::vector<int> unique(const std::vector<int>& u) {
  std::vector<int> v;
  v.reserve(u.size());
  v.push_back(u[0]);
  for (int j = 1; j < u.size(); j++) {
    if (u[j-1] != u[j]) {
      v.push_back(u[j]);
    }
  }
  return v;
}

const char Dna(int i) {
  static const char *dna = "ACGT";
  return dna[i & 0x03];
}

int hamming(const char *a, const char *b) {
  int h = 0;
  while (*a != 0 && *b != 0) {
    if (*a != *b) {
      h++;
    }
    a++;
    b++;
  }
  return h;
}

std::string revcomp(const std::string s) {
  std::string r(s);
  std::transform(s.rbegin(), s.rend(), r.begin(), [](char c) {
      switch(c) {
      case 'A': return 'T';
      case 'C': return 'G';
      case 'G': return 'C';
      case 'T': return 'A';
      default: return 'N';
      }
      return 'N';
    });
  return r;
}

void KmerIndex::Build(const ProgramOptions& opt) {
  // read input
  std::unordered_set<std::string> unique_names;
  int k = opt.k;
  for (auto& fasta : opt.transfasta) {
    std::cerr << "[build] loading fasta file " << fasta
              << std::endl;
  }
  std::cerr << "[build] k-mer length: " << k << std::endl;


  std::vector<std::string> seqs;

  // read fasta file  
  gzFile fp = 0;
  kseq_t *seq;
  int l = 0;
  std::mt19937 gen(42);
  int countNonNucl = 0;
  int countUNuc = 0;
  int polyAcount = 0;

  for (auto& fasta : opt.transfasta) {
    fp = gzopen(fasta.c_str(), "r");
    seq = kseq_init(fp);
    while (true) {
      l = kseq_read(seq);
      if (l <= 0) {
        break;
      }
      seqs.emplace_back(seq->seq.s);
      std::string& str = *seqs.rbegin();
      auto n = str.size();
      for (auto i = 0; i < n; i++) {
        char c = str[i];
        c = ::toupper(c);
        if (c=='U') {
          str[i] = 'T';
          countUNuc++;
        } else if (c !='A' && c != 'C' && c != 'G' && c != 'T') {
          str[i] = Dna(gen()); // replace with pseudorandom string
          countNonNucl++;
        }
      }
      std::transform(str.begin(), str.end(),str.begin(), ::toupper);

      if (str.size() >= 10 && str.substr(str.size()-10,10) == "AAAAAAAAAA") {
        // clip off polyA tail
        //std::cerr << "[index] clipping off polyA tail" << std::endl;
        polyAcount++;
        int j;
        for (j = str.size()-1; j >= 0 && str[j] == 'A'; j--) {}
        str = str.substr(0,j+1);
      }

    
      target_lens_.push_back(seq->seq.l);
      std::string name(seq->name.s);
      size_t p = name.find(' ');
      if (p != std::string::npos) {
        name = name.substr(0,p);
      }

      if (unique_names.find(name) != unique_names.end()) {
        for (int i = 1; ; i++) { // potential bug if you have more than 2^32 repeated names
          std::string new_name = name + "_" + std::to_string(i);
          if (unique_names.find(new_name) == unique_names.end()) {
            name = new_name;
            break;
          }
        }
      }
      unique_names.insert(name);
      target_names_.push_back(name);

    }
    gzclose(fp);
    fp=0;
  }

  if (polyAcount > 0) {
    std::cerr << "[build] warning: clipped off poly-A tail (longer than 10)" << std::endl << "        from " << polyAcount << " target sequences" << std::endl;
  }

  
  if (countNonNucl > 0) {
    std::cerr << "[build] warning: replaced " << countNonNucl << " non-ACGUT characters in the input sequence" << std::endl << "        with pseudorandom nucleotides" << std::endl;
  }
  if (countUNuc > 0) {
    std::cerr << "[build] warning: replaced " << countUNuc << " U characters with Ts" << std::endl;
  }
  
  num_trans = seqs.size();
  
  // for each target, create it's own equivalence class
  for (int i = 0; i < seqs.size(); i++ ) {
    std::vector<int> single(1,i);
    //ecmap.insert({i,single});
    ecmap.push_back(single);
    ecmapinv.insert({single,i});
  }
  
  Buildhashtable(seqs);
  hashtable_aligner();
  write(opt.index);
}

void KmerIndex::hashtable_aligner() {
	std::cerr << "[build] align table size  ... "; std::cerr.flush();
	roundup_4(t_max);
	for(auto& table: *hash_tables) {
		for(auto& entry: *table) {
			while(entry.second.size() < t_max) {
				entry.second.push_back(-1);
			}
		}
		if(table->size() < kmer_max) {
			std::vector<int16_t> tmp(t_max, -1);
			Kmer tmp_kmer; 
			tmp_kmer.set_empty();
			table->insert(std::make_pair(tmp_kmer,-1));
		}
	}
	std::cerr << " done " << std::endl;
}

void KmerIndex::Buildhashtable(const std::vector<std::string>& seqs) {
	std::cerr << "[build] build hash table ... "; std::cerr.flush();
	// gather all k-mers for each group
	// #seqs in a group = ((#seqs)/(#dpu in a pipeline worker))
	uint32_t seq_in_group = seqs.size() / dpu_n;
	uint32_t remain = seqs.size() % dpu_n;
	int i = 0;
	for (int d = 0; d < dpu_n; d++) {
		uint32_t seq_n = (d < remain) ? seq_in_group + 1 : seq_in_group;
		std::map<Kmer, std::vector<int16_t>>* group_hash = new std::map<Kmer, std::vector<int16_t>>();
		for (; i < seq_n; i++) {
			const char *s = seqs[i].c_str();
			KmerIterator kit(s),kit_end;
			for (; kit != kit_end; ++kit) {
				auto found = group_hash->find(kit->first);
				if(found != group_hash->end()) {
					std::vector<int16_t> tmp;
					tmp.push_back(i);
					group_hash->insert(std::make_pair(kit->first, tmp));
				}
				else {
					found->second.push_back(i);
					if(t_max < found->second.size())
						t_max = found->second.size();
				}
			}
		}
		hash_tables->push_back(group_hash);
		if(group_hash->size() > kmer_max)
			kmer_max = group_hash->size();
		delete group_hash;
  }
  std::cerr << " done " << std::endl;
}

void KmerIndex::write(const std::string& index_out, bool writeKmerTable) {
	std::ofstream out;
	out.open(index_out, std::ios::out | std::ios::binary);

	if (!out.is_open()) {
		// TODO: better handling
		std::cerr << "Error: index output file could not be opened!";
		exit(1);
	}

	// write k
	out.write((char *)&k, sizeof(k));

	// write hash table size
	out.write((char *)&t_max, sizeof(t_max));
	out.write((char *)&kmer_max, sizeof(kmer_max));

	// write number of dpus
	out.write((char *)&dpu_n, sizeof(dpu_n));

	// build & write hash table and  
	for(auto& table: *hash_tables) {
		// bulid 
		for(auto& entry: *table) {
			kmap.insert(std::make_pair(entry.first, entry.second));
		}
		// write 
		for (auto& kv : kmap) {
			// write key : kmer
			uint64_t key = kv.first.get_kmer();
			out.write((char *)&key, sizeof(key));
			// write value : transcripts
			for(auto& t : kv.second) {
				out.write((char *)&t, sizeof(t));
			}
		}
	}
	out.flush();
	out.close();
}

bool KmerIndex::fwStep(Kmer km, Kmer& end) const {
  int j = -1;
  int fw_count = 0;
  for (int i = 0; i < 4; i++) {
    Kmer fw_rep = end.forwardBase(Dna(i)).rep();
    auto search = kmap.find(fw_rep);
    if (search != kmap.end()) {
      j = i;
      ++fw_count;
      if (fw_count > 1) {
        return false;
      }
    }
  }

  if (fw_count != 1) {
    return false;
  }

  Kmer fw = end.forwardBase(Dna(j));

  int bw_count = 0;
  for (int i = 0; i < 4; i++) {
    Kmer bw_rep = fw.backwardBase(Dna(i)).rep();
    if (kmap.find(bw_rep) != kmap.end()) {
      ++bw_count;
      if (bw_count > 1) {
        return false;
      }
    }
  }

  if (bw_count != 1) {
    return false;
  } else {
    if (fw != km) {
      end = fw;
      return true;
    } else {
      return false;
    }
  }

}

void KmerIndex::load(ProgramOptions& opt) {
	std::string& index_in = opt.index;
	std::ifstream in;


	in.open(index_in, std::ios::in | std::ios::binary);

	if (!in.is_open()) {
		// TODO: better handling
		std::cerr << "Error: index input file could not be opened!";
		exit(1);
	}

	// read k
	in.read((char *)&k, sizeof(k));

	// read size of hash table
	in.read((char *)&t_max, sizeof(t_max));
	in.read((char *)&kmer_max, sizeof(kmer_max));

	// read number of dpus
	in.read((char *)&dpu_n, sizeof(dpu_n));

	for (int i = 0; i < dpu_n; i++) {
		table_buf.push_back(std::vector<uint64_t>());
		// read 
		for (int j = 0; j < kmer_max; j++) {
			// key 
			uint64_t key;
			in.read((char *)&key, sizeof(key));
			// to transfer buffer
			table_buf[i].push_back(key);
			// value : transcripts
			uint64_t values = 0UL;
			int tt = 0;
			for(int t = 0; t < t_max; t++) {
				int16_t tran;
				in.read((char *)&tran, sizeof(tran));
				if(tt < 4) {
					values |= ((uint64_t)tran);
					if(tt < 3)
						values << 16;
					tt ++;
				}
				if(tt == 4) {
					// to transfer buffer
					table_buf[i].push_back(values);
					tt = 0;
					values = 0UL;
				}
			}
		}
		// TODO : transfer to DPU
		kmer_max_buf = std::vector<int32_t>(dpu_n, kmer_max);
		t_max_buf = std::vector<int32_t>(dpu_n, t_max);
		k_buf = std::vector<int32_t>(dpu_n, k);
	}
}


// std::pair<int,bool> KmerIndex::findPosition(int tr, Kmer km, int p) const {
//   auto it = kmap.find(km.rep());
//   if (it != kmap.end()) {
//     KmerEntry val = it->second;
//     return findPosition(tr, km, val, p);
//   } else {
//     return {-1,true};
//   }
// }

// //use:  (pos,sense) = index.findPosition(tr,km,val,p)
// //pre:  index.kmap[km] == val,
// //      km is the p-th k-mer of a read
// //      val.contig maps to tr
// //post: km is found in position pos (1-based) on the sense/!sense strand of tr
// std::pair<int,bool> KmerIndex::findPosition(int tr, Kmer km, KmerEntry val, int p) const {
//   bool fw = (km == km.rep());
//   bool csense = (fw == val.isFw());

//   int trpos = -1;
//   bool trsense = true;
//   if (val.contig < 0) {
//     return {-1, true};
//   }
//   const Contig &c = dbGraph.contigs[val.contig];
//   for (auto x : c.transcripts) {
//     if (x.trid == tr) {
//       trpos = x.pos;
//       trsense = x.sense;
//       break;
//     }
//   }

//   if (trpos == -1) {
//     return {-1,true};
//   }


//   if (trsense) {
//     if (csense) {
//       return {trpos + val.getPos() - p + 1, csense}; // 1-based, case I
//     } else {
//       return {trpos + val.getPos() + k + p, csense}; // 1-based, case III
//     }
//   } else {
//     if (csense) {
//       return {trpos + (c.length - val.getPos() -1) + k + p, !csense};  // 1-based, case IV
//     } else {
//       return {trpos + (c.length - val.getPos())  - p, !csense}; // 1-based, case II
//     }
//   }
// }

// /*
// // use:  r = matchEnd(s,l,v,p)
// // pre:  v is initialized, p>=0
// // post: v contains all equiv classes for the k-mer in s, as
// //       well as the best match for s[p:]
// //       if r is false then no match was found
// bool KmerIndex::matchEnd(const char *s, int l, std::vector<std::pair<int,int>> &v, int maxPos) const {
//   // kmer-iterator checks for N's and out of bounds
//   KmerIterator kit(s+maxPos), kit_end;
//   if (kit != kit_end && kit->second == 0) {
//     Kmer last = kit->first;
//     auto search = kmap.find(last.rep());

//     if (search == kmap.end()) {
//       return false; // shouldn't happen
//     }

//     KmerEntry val = search->second;
//     bool forward = (kit->first == search->first);
//     int dist = val.getDist(forward);
//     int pos = maxPos + dist + 1; // move 1 past the end of the contig
    
//     const char* sLeft = s + pos + (k-1);
//     int sizeleft = l -  (pos + k-1); // characters left in sleft
//     if (sizeleft <= 0) {
//       return false; // nothing left to match
//     }

//     // figure out end k-mer
//     const Contig& root = dbGraph.contigs[val.contig];
//     Kmer end; // last k-mer
//     bool readFw = (forward == val.isFw());
//     if (readFw) {
//       end = Kmer(root.seq.c_str() + root.length-1);
//     } else {
//       end = Kmer(root.seq.c_str()).twin();
//     }

    
//     int bestContig = -1;
//     int bestDist = sizeleft;
//     int numBest = 0;
//     for (int i = 0; i < 4; i++) {
//       Kmer x = end.forwardBase(Dna(i));
//       Kmer xr = x.rep();
//       auto searchx = kmap.find(xr);
//       if (searchx != kmap.end()) {
//         KmerEntry valx = searchx->second;
//         const Contig& branch = dbGraph.contigs[valx.contig];
//         if (branch.length < sizeleft) {
//           return false; // haven't implemented graph walks yet
//         }
//         std::string cs;
//         bool contigFw = (x==xr) && valx.isFw();
//         // todo: get rid of this string copy
//         if (valx.getPos() == 0 && contigFw) {
//           cs = branch.seq.substr(k-1);
//         } else if (valx.getPos() == branch.length-1 && !contigFw) {
//           cs = revcomp(branch.seq).substr(k-1);
//         }
//         int hdist = hamming(sLeft,cs.c_str());
//         if (bestDist >= hdist) {
//           numBest++;
//           if (bestDist > hdist) {
//             bestDist = hdist;
//             numBest = 1;
//           }
//           bestContig = valx.contig;
//         }
//       }
//     }

//     if (numBest == 1 && bestDist < 10) {
//       v.push_back({dbGraph.ecs[bestContig], l-k});
//       return true;
//     } else {
//       return false;
//     }
//   } else {
//     return false;
//   }
// }
// */


// // use:  res = intersect(ec,v)
// // pre:  ec is in ecmap, v is a vector of valid targets
// //       v is sorted in increasing order
// // post: res contains the intersection  of ecmap[ec] and v sorted increasing
// //       res is empty if ec is not in ecma
// std::vector<int> KmerIndex::intersect(int ec, const std::vector<int>& v) const {
//   std::vector<int> res;
//   //auto search = ecmap.find(ec);
//   if (ec < ecmap.size()) {
//     //if (search != ecmap.end()) {
//     //auto& u = search->second;
//     auto& u = ecmap[ec];
//     res.reserve(v.size());

//     auto a = u.begin();
//     auto b = v.begin();

//     while (a != u.end() && b != v.end()) {
//       if (*a < *b) {
//         ++a;
//       } else if (*b < *a) {
//         ++b;
//       } else {
//         // match
//         res.push_back(*a);
//         ++a;
//         ++b;
//       }
//     }
//   }
//   return res;
// }


// void KmerIndex::loadTranscriptSequences() const {
//   if (target_seqs_loaded) {
//     return;
//   }


  
//   std::vector<std::vector<std::pair<int, ContigToTranscript>>> trans_contigs(num_trans);
//   for (auto &c : dbGraph.contigs) {
//     for (auto &ct : c.transcripts) {
//       trans_contigs[ct.trid].push_back({c.id, ct});
//     }
//   }

//   auto &target_seqs = const_cast<std::vector<std::string>&>(target_seqs_);
  
//   for (int i = 0; i < trans_contigs.size(); i++) {
//     auto &v = trans_contigs[i];
//     std::sort(v.begin(), v.end(), [](std::pair<int,ContigToTranscript> a, std::pair<int,ContigToTranscript> b) {
//         return a.second.pos < b.second.pos;
//       });

//     std::string seq;
//     seq.reserve(target_lens_[i]);

//     for (auto &pct : v) {
//       auto ct = pct.second;
//       int start = (ct.pos==0) ? 0 : k-1;
//       const auto& contig = dbGraph.contigs[pct.first];
//       if (ct.sense) {
//         seq.append(contig.seq.substr(start));
//       } else {
//         seq.append(revcomp(contig.seq).substr(start));
//       }
//     }
//     target_seqs.push_back(seq);
//   }

//   bool &t = const_cast<bool&>(target_seqs_loaded);
//   t = true;//target_seqs_loaded = true;
//   return;
// }

// void KmerIndex::clear() {
//   kmap.clear_table();
//   ecmap.resize(0);
//   dbGraph.ecs.resize(0);
//   dbGraph.contigs.resize(0);
//   {
//     std::unordered_map<std::vector<int>, int, SortedVectorHasher> empty;
//     std::swap(ecmapinv, empty);
//   }
  
//   target_lens_.resize(0);
//   target_names_.resize(0);
//   target_seqs_.resize(0);
// }

// void KmerIndex::writePseudoBamHeader(std::ostream &o) const {
//   // write out header
//   o << "@HD\tVN:1.0\n";
//   for (int i = 0; i < num_trans; i++) {
//     o << "@SQ\tSN:" << target_names_[i] << "\tLN:" << target_lens_[i] << "\n";
//   }
//   o << "@PG\tID:kallisto\tPN:kallisto\tVN:"<< KALLISTO_VERSION << "\n";
//   o.flush();
// }
