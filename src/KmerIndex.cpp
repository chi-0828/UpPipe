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
  write(opt.index);
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

  // 1. write version
  out.write((char *)&INDEX_VERSION, sizeof(INDEX_VERSION));

  // 2. write k
  out.write((char *)&k, sizeof(k));

  // 3. write number of targets
  out.write((char *)&num_trans, sizeof(num_trans));

  // 4. write out target lengths
  for (int tlen : target_lens_) {
    out.write((char *)&tlen, sizeof(tlen));
  }

  size_t kmap_size = kmap.size();

  if (writeKmerTable) {
    // 5. write number of k-mers in map
    out.write((char *)&kmap_size, sizeof(kmap_size));

    // 6. write kmer->ec values
    for (auto& kv : kmap) {
      out.write((char *)&kv.first, sizeof(kv.first));
      out.write((char *)&kv.second, sizeof(kv.second));
    }
  } else {
    // 5. write fake k-mer size
    kmap_size = 0;
    out.write((char *)&kmap_size, sizeof(kmap_size));

    // 6. write none of the kmer->ec values
  }
  // 7. write number of equivalence classes
  size_t tmp_size;
  tmp_size = ecmap.size();
  out.write((char *)&tmp_size, sizeof(tmp_size));

  // 8. write out each equiv class
  //  for (auto& kv : ecmap) {
  for (int ec = 0; ec < ecmap.size(); ec++) {
    out.write((char *)&ec, sizeof(ec));
    auto& v = ecmap[ec];
    // 8.1 write out the size of equiv class
    tmp_size = v.size();
    out.write((char *)&tmp_size, sizeof(tmp_size));
    // 8.2 write each member
    for (auto& val: v) {
      out.write((char *)&val, sizeof(val));
    }
  }

  // 9. Write out target ids
  // XXX: num_trans should equal to target_names_.size(), so don't need
  // to write out again.
  assert(num_trans == target_names_.size());
  for (auto& tid : target_names_) {
    // 9.1 write out how many bytes
    // XXX: Note: this doesn't actually encore the max targ id size.
    // might cause problems in the future
    // tmp_size = tid.size();
    tmp_size = strlen(tid.c_str());
    out.write((char *)&tmp_size, sizeof(tmp_size));

    // 9.2 write out the actual string
    out.write(tid.c_str(), tmp_size);
  }

  // 10. write out contigs
  if (writeKmerTable) {
    assert(dbGraph.contigs.size() == dbGraph.ecs.size());
    tmp_size = dbGraph.contigs.size();
    out.write((char*)&tmp_size, sizeof(tmp_size));
    for (auto& contig : dbGraph.contigs) {
      out.write((char*)&contig.id, sizeof(contig.id));
      out.write((char*)&contig.length, sizeof(contig.length));
      tmp_size = strlen(contig.seq.c_str());
      out.write((char*)&tmp_size, sizeof(tmp_size));
      out.write(contig.seq.c_str(), tmp_size);

      // 10.1 write out transcript info
      tmp_size = contig.transcripts.size();
      out.write((char*)&tmp_size, sizeof(tmp_size));
      for (auto& info : contig.transcripts) {
        out.write((char*)&info.trid, sizeof(info.trid));
        out.write((char*)&info.pos, sizeof(info.pos));
        out.write((char*)&info.sense, sizeof(info.sense));
      }
    }
    
    // 11. write out ecs info
    for (auto ec : dbGraph.ecs) {
      out.write((char*)&ec, sizeof(ec));
    }


  } else {
    // write empty dBG
    tmp_size = 0;
    out.write((char*)&tmp_size, sizeof(tmp_size));
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

void KmerIndex::load(ProgramOptions& opt, bool loadKmerTable) {
  if (opt.index.empty() && !loadKmerTable) {
    // Make an index from transcript and EC files
    loadTranscriptsFromFile(opt);
    loadECsFromFile(opt);
    return;
  }

  std::string& index_in = opt.index;
  std::ifstream in;


  in.open(index_in, std::ios::in | std::ios::binary);

  if (!in.is_open()) {
    // TODO: better handling
    std::cerr << "Error: index input file could not be opened!";
    exit(1);
  }

  // 1. read version
  size_t header_version = 0;
  in.read((char *)&header_version, sizeof(header_version));

  if (header_version != INDEX_VERSION) {
    std::cerr << "Error: incompatible indices. Found version " << header_version << ", expected version " << INDEX_VERSION << std::endl
              << "Rerun with index to regenerate";
    exit(1);
  }

  // 2. read k
  in.read((char *)&k, sizeof(k));
  if (Kmer::k == 0) {
    //std::cerr << "[index] no k has been set, setting k = " << k << std::endl;
    Kmer::set_k(k);
    opt.k = k;
  } else if (Kmer::k == k) {
    //std::cerr << "[index] Kmer::k has been set and matches" << k << std::endl;
    opt.k = k;
  } else {
    std::cerr << "Error: Kmer::k was already set to = " << Kmer::k << std::endl
              << "       conflicts with value of k  = " << k << std::endl;
    exit(1);
  }

  // 3. read in number of targets
  in.read((char *)&num_trans, sizeof(num_trans));

  // 4. read in length of targets
  target_lens_.clear();
  target_lens_.reserve(num_trans);

  for (int i = 0; i < num_trans; i++) {
    int tlen;
    in.read((char *)&tlen, sizeof(tlen));
    target_lens_.push_back(tlen);
  }

  // 5. read number of k-mers
  size_t kmap_size;
  in.read((char *)&kmap_size, sizeof(kmap_size));

  std::cerr << "[index] k-mer length: " << k << std::endl;
  std::cerr << "[index] number of targets: " << pretty_num(num_trans)
    << std::endl;
  std::cerr << "[index] number of k-mers: " << pretty_num(kmap_size)
    << std::endl;

  kmap.clear();
  if (loadKmerTable) {
    kmap.reserve(kmap_size,true);
  }

  // 6. read kmer->ec values
  Kmer tmp_kmer;
  KmerEntry tmp_val;
  for (size_t i = 0; i < kmap_size; ++i) {
    in.read((char *)&tmp_kmer, sizeof(tmp_kmer));
    in.read((char *)&tmp_val, sizeof(tmp_val));

    if (loadKmerTable) {
      kmap.insert({tmp_kmer, tmp_val});
    }
  }

  // 7. read number of equivalence classes
  size_t ecmap_size;
  in.read((char *)&ecmap_size, sizeof(ecmap_size));

  std::cerr << "[index] number of equivalence classes: "
    << pretty_num(ecmap_size) << std::endl;
  ecmap.resize(ecmap_size);
  int tmp_id;
  int tmp_ecval;
  size_t vec_size;
  // 8. read each equiv class
  for (size_t ec = 0; ec < ecmap_size; ++ec) {
    in.read((char *)&tmp_id, sizeof(tmp_id));

    // 8.1 read size of equiv class
    in.read((char *)&vec_size, sizeof(vec_size));

    // 8.2 read each member
    std::vector<int> tmp_vec;
    tmp_vec.reserve(vec_size);
    for (size_t j = 0; j < vec_size; ++j ) {
      in.read((char *)&tmp_ecval, sizeof(tmp_ecval));
      tmp_vec.push_back(tmp_ecval);
    }
    //ecmap.insert({tmp_id, tmp_vec});
    ecmap[tmp_id] = tmp_vec;
    ecmapinv.insert({tmp_vec, tmp_id});
  }

  // 9. read in target ids
  target_names_.clear();
  target_names_.reserve(num_trans);

  size_t tmp_size;
  size_t bufsz = 1024;
  char *buffer = new char[bufsz];
  for (auto i = 0; i < num_trans; ++i) {
    // 9.1 read in the size
    in.read((char *)&tmp_size, sizeof(tmp_size));

    if (tmp_size +1 > bufsz) {
      delete[] buffer;
      bufsz = 2*(tmp_size+1);
      buffer = new char[bufsz];
    }
    
    // clear the buffer 
    memset(buffer,0,bufsz);
    // 9.2 read in the character string
    in.read(buffer, tmp_size);

    /* std::string tmp_targ_id( buffer ); */
    target_names_.push_back(std::string( buffer ));
  }


  // 10. read contigs
  size_t contig_size;
  in.read((char *)&contig_size, sizeof(contig_size));
  dbGraph.contigs.clear();
  dbGraph.contigs.reserve(contig_size);
  for (auto i = 0; i < contig_size; i++) {
    Contig c;
    in.read((char *)&c.id, sizeof(c.id));
    in.read((char *)&c.length, sizeof(c.length));
    in.read((char *)&tmp_size, sizeof(tmp_size));

    if (tmp_size + 1 > bufsz) {
      delete[] buffer;
      bufsz = 2*(tmp_size+1);
      buffer = new char[bufsz];
    }

    memset(buffer,0,bufsz);
    in.read(buffer, tmp_size);
    c.seq = std::string(buffer); // copy
    
    // 10.1 read transcript info
    in.read((char*)&tmp_size, sizeof(tmp_size));
    c.transcripts.clear();
    c.transcripts.reserve(tmp_size);

    for (auto j = 0; j < tmp_size; j++) {
      ContigToTranscript info;
      in.read((char*)&info.trid, sizeof(info.trid));
      in.read((char*)&info.pos, sizeof(info.pos));
      in.read((char*)&info.sense, sizeof(info.sense));
      c.transcripts.push_back(info);
    }

    dbGraph.contigs.push_back(c);
  }

  // 11. read ecs info
  dbGraph.ecs.clear();
  dbGraph.ecs.reserve(contig_size);
  int tmp_ec;
  for (auto i = 0; i < contig_size; i++) {
    in.read((char *)&tmp_ec, sizeof(tmp_ec));
    dbGraph.ecs.push_back(tmp_ec);
  }

  // delete the buffer
  delete[] buffer;
  buffer=nullptr;
  
  in.close();
  
  if (!opt.ecFile.empty()) {
    loadECsFromFile(opt);
  }


  // process for dpu
  if(kmap.size_ > MAX_table_n){
    std::vector<int64_t> table_int;
    std::vector<uint64_t> table_kmer;
    //table_kmer_buf.reserve(opt.dpu);
    //table_int_buf.reserve(opt.dpu);
    //round_buf.reserve(opt.dpu);
    int round = (kmap.size_/MAX_table_n);
    //std::cerr << "round " << round << "\n";
    int head = 0;
    for(int r = 0; r < round; r++){
      size_vec.push_back(std::vector<size_t>(1, MAX_table_n));

      table_int.clear();
      table_kmer.clear();
      for(int i = 0; i<MAX_table_n; i++){
        auto table_tmp = kmap.table;
        table_tmp = table_tmp + head + i;
        table_int.push_back(table_tmp->second.contig);
        table_kmer.push_back(table_tmp->first.get_kmer()[0]);
      }
      table_int_vec.push_back(table_int);
      table_kmer_vec.push_back(table_kmer);
      head += MAX_table_n;
    }
    int remain = (kmap.size_%MAX_table_n);
    //std::cerr << "remain " << remain << "\n";
    //getchar();
    if( remain > 0){
      round ++;
      size_vec.push_back(std::vector<size_t>(1, remain));

      table_int.clear();
      table_kmer.clear();
      for(int i = 0; i<MAX_table_n; i++){
        auto table_tmp = kmap.table;
        if(i< remain){
          table_tmp = table_tmp + head + i;
          table_int.push_back(table_tmp->second.contig);
          table_kmer.push_back(table_tmp->first.get_kmer()[0]);
        }
        else{
          auto table_tmp = kmap.empty;
          table_int.push_back(table_tmp.second.contig);
          table_kmer.push_back(table_tmp.first.get_kmer()[0]);
        }
      }
      table_int_vec.push_back(table_int);
      table_kmer_vec.push_back(table_kmer);
      head += remain;
    }

    // check the table 
    assert(head == kmap.size_);
    // prepare table transfer
    int dpu_round = 0;
    for(int d = 0; d < opt.dpu_n; d++){
      table_kmer_buf.push_back(std::vector<size_t>());
      table_int_buf.push_back(std::vector<int64_t>());
      round_buf.push_back(std::vector<int32_t>());
      round_buf[d].push_back(dpu_round);
      table_int_buf[d].insert(table_int_buf[d].end(), table_int_vec[dpu_round].begin(), table_int_vec[dpu_round].end());
      table_kmer_buf[d].insert(table_kmer_buf[d].end(), table_kmer_vec[dpu_round].begin(), table_kmer_vec[dpu_round].end());
      dpu_round++ ;
      if(dpu_round == round){
        dpu_round = 0;
      }
    }
  }
}

void KmerIndex::loadECsFromFile(const ProgramOptions& opt) {
  ecmap.clear();
  ecmapinv.clear();
  int i = 0;
  std::ifstream in((opt.ecFile));
  if (in.is_open()) {
    std::string line;
    while (getline(in, line)) {
      std::stringstream ss(line);
      int ec;
      std::string transcripts;
      ss >> ec >> transcripts;
      if (i != ec) {
        std::cerr << "Error: equivalence class file has a misplaced equivalence class."
                  << " Found " << ec << ", expected " << i << std::endl;
        exit(1);
      }
      std::vector<int> tmp_vec;
      std::stringstream ss2(transcripts);
      while(ss2.good()) {
        std::string tmp_ecval;
        getline(ss2, tmp_ecval, ',');
        int tmp_ecval_num = std::atoi(tmp_ecval.c_str());
        if (tmp_ecval_num < 0 || tmp_ecval_num >= num_trans) {
          std::cerr << "Error: equivalence class file has invalid value: " 
                    << tmp_ecval << " in " << transcripts << std::endl;
          exit(1);
        }
        tmp_vec.push_back(tmp_ecval_num);
      }
      ecmap.push_back(tmp_vec); // copy
      ecmapinv.insert({std::move(tmp_vec), i}); // move
      i++;
    }
  } else {
    std::cerr << "Error: could not open file " << opt.ecFile << std::endl;
    exit(1);
  }
  std::cerr << "[index] number of equivalence classes loaded from file: "
            << pretty_num(ecmap.size()) << std::endl;
}

void KmerIndex::loadTranscriptsFromFile(const ProgramOptions& opt) {
  target_names_.clear();
  int i = 0;
  std::ifstream in((opt.transcriptsFile));
  if (in.is_open()) {
    std::string txp;
    while (in >> txp) {
      target_names_.push_back(txp);
      i++;
    }
  } else {
    std::cerr << "Error: could not open file " << opt.transcriptsFile << std::endl;
    exit(1);
  }
  num_trans = i;
  target_lens_.assign(num_trans, 0);
  std::cerr << "[index] number of targets loaded from file: "
            << pretty_num(num_trans) << std::endl;
}

int KmerIndex::mapPair(const char *s1, int l1, const char *s2, int l2, int ec) const {
  bool d1 = true;
  bool d2 = true;
  int p1 = -1;
  int p2 = -1;
  int c1 = -1;
  int c2 = -1;


  KmerIterator kit1(s1), kit_end;

  bool found1 = false;
  for (; kit1 != kit_end; ++kit1) {
    Kmer x = kit1->first;
    Kmer xr = x.rep();
    auto search = kmap.find(xr);
    bool forward = (x==xr);

    if (search != kmap.end()) {
      found1 = true;
      KmerEntry val = search->second;
      c1 = val.contig;
      if (forward == val.isFw()) {
        p1 = val.getPos() - kit1->second;
        d1 = true;
      } else {
        p1 = val.getPos() + k + kit1->second;
        d1 = false;
      }
      break;
    }
  }

  if (!found1) {
    return -1;
  }

  

  KmerIterator kit2(s2);
  bool found2 = false;

  for (; kit2 != kit_end; ++kit2) {
    Kmer x = kit2->first;
    Kmer xr = x.rep();
    auto search = kmap.find(xr);
    bool forward = (x==xr);

    if (search != kmap.end()) {
      found2 = true;
      KmerEntry val = search->second;
      c2 = val.contig;
      if (forward== val.isFw()) {
        p2 = val.getPos() - kit2->second;
        d2 = true;
      } else {
        p2 = val.getPos() + k + kit2->second;
        d2 = false;
      }
      break;
    }
  }

  if (!found2) {
    return -1;
  }

  if (c1 != c2) {
    return -1;
  }

  if ((d1 && d2) || (!d1 && !d2)) {
    //std::cerr << "Reads map to same strand " << s1 << "\t" << s2 << std::endl;
    return -1;
  }

  if (p1>p2) {
    return p1-p2;
  } else {
    return p2-p1;
  }

}

// use:  match(s,l,v)
// pre:  v is initialized
// post: v contains all equiv classes for the k-mers in s
void KmerIndex::match(const char *s, int l, std::vector<std::pair<KmerEntry, int>>& v) const {
  KmerIterator kit(s), kit_end;
  bool backOff = false;
  int nextPos = 0; // nextPosition to check
  for (int i = 0;  kit != kit_end; ++i,++kit) {
    // need to check it
    auto search = kmap.find(kit->first.rep());
    int pos = kit->second;

    if (search != kmap.end()) {

      KmerEntry val = search->second;
      
      //std::cerr << "1148 push: " << val.contig << "(kmer: " << kit->first.toString() << ")\n";
      v.push_back({val, kit->second});

      // see if we can skip ahead
      // bring thisback later
      bool forward = (kit->first == search->first);
      int dist = val.getDist(forward);


      //const int lastbp = 10;
      if (dist >= 2) {
        // where should we jump to?
        int nextPos = pos+dist; // default jump

        if (pos + dist >= l-k) {
          // if we can jump beyond the read, check the end
          nextPos = l-k;
        }

        // check next position
        KmerIterator kit2(kit);
        kit2.jumpTo(nextPos);
        if (kit2 != kit_end) {
          Kmer rep2 = kit2->first.rep();
          auto search2 = kmap.find(rep2);
          bool found2 = false;
          int  found2pos = pos+dist;
          if (search2 == kmap.end()) {
            found2=true;
            found2pos = pos;
          } else if (val.contig == search2->second.contig) {
            found2=true;
            found2pos = pos+dist;
          }
          if (found2) {
            // great, a match (or nothing) see if we can move the k-mer forward
            if (found2pos >= l-k) {
              //std::cerr << "1185 push: " << val.contig << "(kmer: " << kit2->first.toString() << ")\n";
              v.push_back({val, l-k}); // push back a fake position
              break; //
            } else {
              //std::cerr << "1189 push: " << val.contig << "(kmer: " << kit2->first.toString() << ")\n";
              v.push_back({val, found2pos});
              kit = kit2; // move iterator to this new position
            }
          } else {
            // this is weird, let's try the middle k-mer
            bool foundMiddle = false;
            if (dist > 4) {
              int middlePos = (pos + nextPos)/2;
              int middleContig = -1;
              int found3pos = pos+dist;
              KmerIterator kit3(kit);
              kit3.jumpTo(middlePos);
              KmerEntry val3;
              if (kit3 != kit_end) {
                Kmer rep3 = kit3->first.rep();
                auto search3 = kmap.find(rep3);
                if (search3 != kmap.end()) {
                  middleContig = search3->second.contig;
                  if (middleContig == val.contig) {
                    foundMiddle = true;
                    found3pos = middlePos;
                  } else if (middleContig == search2->second.contig) {
                    foundMiddle = true;
                    found3pos = pos+dist;
                  }
                }


                if (foundMiddle) {
                  //std::cerr << "1219 push: " << search3->second.contig << "(kmer: " << kit3->first.toString() << ")\n";
                  v.push_back({search3->second, found3pos});
                  if (nextPos >= l-k) {
                    break;
                  } else {
                    kit = kit2; 
                  }
                }
              }
            }


            if (!foundMiddle) {
              ++kit;
              backOff = true;
              goto donejumping; // sue me Dijkstra!
            }
          }
        } else {
          // the sequence is messed up at this point, let's just take the match
          //v.push_back({dbGraph.ecs[val.contig], l-k});
          break;
        }
      }
    }

donejumping:

    if (backOff) {
      // backup plan, let's play it safe and search incrementally for the rest, until nextStop
      for (int j = 0; kit != kit_end; ++kit,++j) {
        if (j==1) {
          j=0;
        }
        if (j==0) {
          // need to check it
          Kmer rep = kit->first.rep();
          auto search = kmap.find(rep);
          if (search != kmap.end()) {
            // if k-mer found
            //std::cerr << "1259 push: " << search->second.contig << "(kmer: " << kit->first.toString() << ")\n";
            v.push_back({search->second, kit->second}); // add equivalence class, and position
          }
        }

        if (kit->second >= nextPos) {
          backOff = false;
          break; // break out of backoff for loop
        }
      }
    }
  }
}

std::pair<int,bool> KmerIndex::findPosition(int tr, Kmer km, int p) const {
  auto it = kmap.find(km.rep());
  if (it != kmap.end()) {
    KmerEntry val = it->second;
    return findPosition(tr, km, val, p);
  } else {
    return {-1,true};
  }
}

//use:  (pos,sense) = index.findPosition(tr,km,val,p)
//pre:  index.kmap[km] == val,
//      km is the p-th k-mer of a read
//      val.contig maps to tr
//post: km is found in position pos (1-based) on the sense/!sense strand of tr
std::pair<int,bool> KmerIndex::findPosition(int tr, Kmer km, KmerEntry val, int p) const {
  bool fw = (km == km.rep());
  bool csense = (fw == val.isFw());

  int trpos = -1;
  bool trsense = true;
  if (val.contig < 0) {
    return {-1, true};
  }
  const Contig &c = dbGraph.contigs[val.contig];
  for (auto x : c.transcripts) {
    if (x.trid == tr) {
      trpos = x.pos;
      trsense = x.sense;
      break;
    }
  }

  if (trpos == -1) {
    return {-1,true};
  }


  if (trsense) {
    if (csense) {
      return {trpos + val.getPos() - p + 1, csense}; // 1-based, case I
    } else {
      return {trpos + val.getPos() + k + p, csense}; // 1-based, case III
    }
  } else {
    if (csense) {
      return {trpos + (c.length - val.getPos() -1) + k + p, !csense};  // 1-based, case IV
    } else {
      return {trpos + (c.length - val.getPos())  - p, !csense}; // 1-based, case II
    }
  }
}

/*
// use:  r = matchEnd(s,l,v,p)
// pre:  v is initialized, p>=0
// post: v contains all equiv classes for the k-mer in s, as
//       well as the best match for s[p:]
//       if r is false then no match was found
bool KmerIndex::matchEnd(const char *s, int l, std::vector<std::pair<int,int>> &v, int maxPos) const {
  // kmer-iterator checks for N's and out of bounds
  KmerIterator kit(s+maxPos), kit_end;
  if (kit != kit_end && kit->second == 0) {
    Kmer last = kit->first;
    auto search = kmap.find(last.rep());

    if (search == kmap.end()) {
      return false; // shouldn't happen
    }

    KmerEntry val = search->second;
    bool forward = (kit->first == search->first);
    int dist = val.getDist(forward);
    int pos = maxPos + dist + 1; // move 1 past the end of the contig
    
    const char* sLeft = s + pos + (k-1);
    int sizeleft = l -  (pos + k-1); // characters left in sleft
    if (sizeleft <= 0) {
      return false; // nothing left to match
    }

    // figure out end k-mer
    const Contig& root = dbGraph.contigs[val.contig];
    Kmer end; // last k-mer
    bool readFw = (forward == val.isFw());
    if (readFw) {
      end = Kmer(root.seq.c_str() + root.length-1);
    } else {
      end = Kmer(root.seq.c_str()).twin();
    }

    
    int bestContig = -1;
    int bestDist = sizeleft;
    int numBest = 0;
    for (int i = 0; i < 4; i++) {
      Kmer x = end.forwardBase(Dna(i));
      Kmer xr = x.rep();
      auto searchx = kmap.find(xr);
      if (searchx != kmap.end()) {
        KmerEntry valx = searchx->second;
        const Contig& branch = dbGraph.contigs[valx.contig];
        if (branch.length < sizeleft) {
          return false; // haven't implemented graph walks yet
        }
        std::string cs;
        bool contigFw = (x==xr) && valx.isFw();
        // todo: get rid of this string copy
        if (valx.getPos() == 0 && contigFw) {
          cs = branch.seq.substr(k-1);
        } else if (valx.getPos() == branch.length-1 && !contigFw) {
          cs = revcomp(branch.seq).substr(k-1);
        }
        int hdist = hamming(sLeft,cs.c_str());
        if (bestDist >= hdist) {
          numBest++;
          if (bestDist > hdist) {
            bestDist = hdist;
            numBest = 1;
          }
          bestContig = valx.contig;
        }
      }
    }

    if (numBest == 1 && bestDist < 10) {
      v.push_back({dbGraph.ecs[bestContig], l-k});
      return true;
    } else {
      return false;
    }
  } else {
    return false;
  }
}
*/


// use:  res = intersect(ec,v)
// pre:  ec is in ecmap, v is a vector of valid targets
//       v is sorted in increasing order
// post: res contains the intersection  of ecmap[ec] and v sorted increasing
//       res is empty if ec is not in ecma
std::vector<int> KmerIndex::intersect(int ec, const std::vector<int>& v) const {
  std::vector<int> res;
  //auto search = ecmap.find(ec);
  if (ec < ecmap.size()) {
    //if (search != ecmap.end()) {
    //auto& u = search->second;
    auto& u = ecmap[ec];
    res.reserve(v.size());

    auto a = u.begin();
    auto b = v.begin();

    while (a != u.end() && b != v.end()) {
      if (*a < *b) {
        ++a;
      } else if (*b < *a) {
        ++b;
      } else {
        // match
        res.push_back(*a);
        ++a;
        ++b;
      }
    }
  }
  return res;
}


void KmerIndex::loadTranscriptSequences() const {
  if (target_seqs_loaded) {
    return;
  }


  
  std::vector<std::vector<std::pair<int, ContigToTranscript>>> trans_contigs(num_trans);
  for (auto &c : dbGraph.contigs) {
    for (auto &ct : c.transcripts) {
      trans_contigs[ct.trid].push_back({c.id, ct});
    }
  }

  auto &target_seqs = const_cast<std::vector<std::string>&>(target_seqs_);
  
  for (int i = 0; i < trans_contigs.size(); i++) {
    auto &v = trans_contigs[i];
    std::sort(v.begin(), v.end(), [](std::pair<int,ContigToTranscript> a, std::pair<int,ContigToTranscript> b) {
        return a.second.pos < b.second.pos;
      });

    std::string seq;
    seq.reserve(target_lens_[i]);

    for (auto &pct : v) {
      auto ct = pct.second;
      int start = (ct.pos==0) ? 0 : k-1;
      const auto& contig = dbGraph.contigs[pct.first];
      if (ct.sense) {
        seq.append(contig.seq.substr(start));
      } else {
        seq.append(revcomp(contig.seq).substr(start));
      }
    }
    target_seqs.push_back(seq);
  }

  bool &t = const_cast<bool&>(target_seqs_loaded);
  t = true;//target_seqs_loaded = true;
  return;
}

void KmerIndex::clear() {
  kmap.clear_table();
  ecmap.resize(0);
  dbGraph.ecs.resize(0);
  dbGraph.contigs.resize(0);
  {
    std::unordered_map<std::vector<int>, int, SortedVectorHasher> empty;
    std::swap(ecmapinv, empty);
  }
  
  target_lens_.resize(0);
  target_names_.resize(0);
  target_seqs_.resize(0);
}

void KmerIndex::writePseudoBamHeader(std::ostream &o) const {
  // write out header
  o << "@HD\tVN:1.0\n";
  for (int i = 0; i < num_trans; i++) {
    o << "@SQ\tSN:" << target_names_[i] << "\tLN:" << target_lens_[i] << "\n";
  }
  o << "@PG\tID:kallisto\tPN:kallisto\tVN:"<< KALLISTO_VERSION << "\n";
  o.flush();
}
