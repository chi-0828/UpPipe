#include "KmerIndex.h"
#include <algorithm>
#include <random>
#include <sstream>
#include <ctype.h>
#include <zlib.h>
#include <unordered_set>
#include "kseq.h"
#include <chrono>


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
  std::cerr << "[build] loading fasta file " << opt.transfasta << std::endl;
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


  fp = gzopen(opt.transfasta.c_str(), "r");
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
  
//   // for each target, create it's own equivalence class
//   for (int i = 0; i < seqs.size(); i++ ) {
//     std::vector<int> single(1,i);
//     //ecmap.insert({i,single});
//     ecmap.push_back(single);
//     ecmapinv.insert({single,i});
//   }
  
  Buildhashtable(seqs);
  hashtable_aligner();
  write(opt.index);
}

void KmerIndex::hashtable_aligner() {
	std::cerr << "[build] align table size  ... "; std::cerr.flush();
	for(auto& table: table_buf) {
		if(kmer_max < table.size())
			kmer_max = table.size();
	}
	for(auto& table: table_buf) {
		for(;table.size() < kmer_max;) {
			table.push_back(EMPTY_KMER);
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
	int16_t i = 0;
	Kmer::set_k(k);
	first_tid_buf.resize(dpu_n);
	table_buf.resize(dpu_n);
	t_max_buf.resize(dpu_n);
	size_buf.resize(dpu_n);
	for (int d = 0; d < dpu_n; d++) {
		uint32_t seq_n = (d < remain) ? seq_in_group + 1 : seq_in_group;
		t_max = 1;
		kmap.init_table(1024);
		int16_t start = i;
		first_tid_buf[d].push_back((int32_t)start);
		for (; ((i < start + (int16_t)seq_n) && i < seqs.size()) ; i++) {
			const char *s = seqs[i].c_str();
			KmerIterator kit(s), kit_end;
			std::deque<Kmer> minimizer;
			for (; kit != kit_end; ++kit) {
				// minimizer
				minimizer.push_back(kit->first);
				if(minimizer.size() >= W(k)) {
					assert(minimizer.size() == W(k));
					// get the mini one
					auto ret = kmap.insert(get_mini(minimizer), i);
					if(t_max < ret)
						t_max = ret;
					minimizer.pop_front();
				}
			}
		}
		roundup_4(t_max);
		t_max_buf[d].push_back(t_max);
		size_buf[d].push_back(kmap.size_);
		// std::cerr << "t_max " << t_max << " kmap.size_" << kmap.size_ << "\n";
		for (size_t h = 0; h < kmap.size_; h++) {
			table_buf[d].push_back(kmap.table[h].first.tobinary());
			int t_size = 0;
			uint64_t v = 0UL;
			// std::cerr <<  kmap.table[h].first.toString()  << " : ";
			for(auto& t : kmap.table[h].second) {
				// std::cerr <<  t  << " ";
				v |= ((uint64_t)t);
				t_size ++;
				if(t_size % 4 == 0) {
					table_buf[d].push_back(v);
					v = 0UL;
				}
				else {
					v <<= 16;
				}
			}
			for(; t_size < t_max;) {
				v |= ((uint64_t)0xFFFF);
				t_size ++;
				if(t_size % 4 == 0) {
					table_buf[d].push_back(v);
					v = 0UL;
				}
				else {
					v <<= 16;
				}
			}
		}
		// uncomment the following code block can print the hash table
		// bool key = 1;
		// int i = 0;
		// int ee = 0;
		// int h = 0;
		// while(ee < table_buf[d].size() ) {
		// 	if(key)
		// 		h ++;
		// 	uint64_t entry = table_buf[d][ee];
		// 	if(key && entry == EMPTY_KMER) {
		// 		// ee += 1;
		// 		ee += (t_max / 4) ;
		// 	} 
		// 	else {
		// 		if(key) {
		// 			std::cerr << "at h = " << h - 1  << " "; 
		// 			std::bitset<64> binary(entry);
		// 			std::cerr << entry << "(" << toString(entry, k)<< ") | ";
		// 			key = 0;
		// 		}
		// 		else {
		// 			int16_t t1 = (int16_t)((((uint64_t)entry)>>(0))&(0xFFFF));
		// 			int16_t t2 = (int16_t)((((uint64_t)entry)>>(16))&(0xFFFF));
		// 			int16_t t3 = (int16_t)((((uint64_t)entry)>>(32))&(0xFFFF));
		// 			int16_t t4 = (int16_t)((((uint64_t)entry)>>(48))&(0xFFFF));
		// 			std::cerr << std::dec << t4 << " " << t3 << " " << t2 << " " << t1 << " ";
		// 			i += 4;
		// 			if( i == t_max) {
		// 				i = 0;
		// 				key = 1;
		// 				std::cerr << "\n";
		// 				//getchar();
		// 			}
		// 		}
		// 	}
		// 	ee ++;
		// }
  	}
  	std::cerr << " done " << std::endl;
}

void KmerIndex::write(const std::string& index_out, bool writeKmerTable) {
	std::cerr << "[build] write hash table to \"" << index_out << "\" ... "; std::cerr.flush();
	std::ofstream out;
	int buf_size = 1024*1024*2;
	char buffer[buf_size];
  	out.rdbuf()->pubsetbuf(buffer, buf_size);
	out.open(index_out, std::ios::out | std::ios::binary);

	if (!out.is_open()) {
		// TODO: better handling
		std::cerr << "Error: index output file \"" << index_out << "\" could not be opened!\n";
		exit(1);
	}

	// write k
	out.write((char *)&k, sizeof(k));
	//std::cerr << k << "\n";

	// write hash table size
	
	out.write((char *)&kmer_max, sizeof(kmer_max));
	//std::cerr << t_max << "\n";
	//std::cerr << kmer_max << "\n";

	// write number of dpus
	out.write((char *)&dpu_n, sizeof(dpu_n));

	// build & write hash table and  
	int d = 0;
	for(auto& table: table_buf) {
		out.write((char *)&first_tid_buf[d][0], sizeof(first_tid_buf[d][0]));
		out.write((char *)&t_max_buf[d][0], sizeof(t_max_buf[d][0]));
		out.write((char *)&size_buf[d][0], sizeof(size_buf[d][0]));
		for(auto& entry: table) {
			out.write((char *)&entry, sizeof(entry));
		}
		d ++;
	}
	out.flush();
	out.close();
	std::cerr << " done " << std::endl;
}

void KmerIndex::load(ProgramOptions& opt) {
	std::cerr << "[load] loading index file " << opt.index << " ... ";
  	auto t_start = std::chrono::high_resolution_clock::now();
	std::string& index_in = opt.index;
	std::ifstream in;

	int buf_size = 1024*1024*2;
	char buffer[buf_size];
  	in.rdbuf()->pubsetbuf(buffer, buf_size); 
	
	in.open(index_in, std::ios::in | std::ios::binary);

	if (!in.is_open()) {
		// TODO: better handling
		std::cerr << "Error: index input file could not be opened!\n";
		exit(1);
	}
	// read chunk
	int buf[3];
	in.read((char *)&buf, sizeof(int)*3);
	// k
	k = buf[0];
	// size of hash table
	kmer_max = buf[1];
	// number of dpus
	dpu_n = buf[2];
	
  	kmer_max_buf = std::vector<int32_t>(1, kmer_max);
	k_buf = std::vector<int32_t>(1, k);

	table_buf.resize(dpu_n);
	t_max_buf.resize(dpu_n);
	size_buf.resize(dpu_n);
	first_tid_buf.resize(dpu_n);

	for (int i = 0; i < dpu_n; i++) {
		int32_t first_tid;
		in.read((char *)&first_tid, sizeof(first_tid));
		first_tid_buf[i].push_back(first_tid);
		int t_max;
		in.read((char *)&t_max, sizeof(t_max));
		t_max_buf[i].push_back(t_max);
		size_t size_;
		in.read((char *)&size_, sizeof(size_));
		size_buf[i].push_back(size_);
		// read 
		for (int j = 0; j < kmer_max/1024; j++) {
			// entry 
			uint64_t buff[1024];
			in.read((char *)&buff, sizeof(uint64_t)*1024);
			for(int a = 0; a < 1024; a ++)
				table_buf[i].push_back(buff[a]);
		}
		for (int j = 0; j < kmer_max%1024; j++) {
			// entry 
			uint64_t entry;
			in.read((char *)&entry, sizeof(entry));
			// to transfer buffer
			table_buf[i].push_back(entry);
		}
	}
	in.close();
	auto t_end = std::chrono::high_resolution_clock::now();
	double get_read_time = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
	std::cerr << "done in " << get_read_time << " s\n";
}

