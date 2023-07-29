#include "dpu_def.h"
#include "hash.h"
#include <mram.h>
#include <assert.h>
#include "defs.h"
#include "mutex.h"
#include <barrier.h>

// about 1MB read packet (from hsot)
__mram_noinit char reads[PACKET_CAPACITY];

// result (to host)
__mram_noinit dpu_result result[PACKET_SIZE];

// 60 MB hash table (from hsot)
__mram_noinit uint64_t table[MAX_table_n];

// read info (from hsot)
__host int32_t reads_len[PACKET_SIZE];
__host int32_t read_n;
__host size_tt size_;

// result info (to hsot)
__host int matched_size;

// runing info (k, t_max, kmer_max)
__host int k;
__host int kmer_max;
__host int t_max;

//STDOUT_BUFFER_INIT(40960)

int RoundDown(int* a)
{
    return *(a) & (-8);
}

int main(){

	// init
	unsigned int tasklet_id = me();
	int gap = t_max/4 + 1;
	// rseult cahce
	//   int64_t *vint = result_id_tasklet[tasklet_id];
	//   int64_t *vpos= result_pos_tasklet[tasklet_id];
	//   v_16len[tasklet_id] = 0;
	//   int v_matched = 0;
	//   int v_len = v_16len[tasklet_id];

	for(int readid = tasklet_id; readid < PACKET_SIZE; readid+=NR_TASKLETS){
		// init result
		__dma_aligned dpu_result t;
		memset(&t, 0, sizeof(dpu_result));
		// read out the read from mram
		int len = READ_LEN;
		__dma_aligned char read_cache[160];
		int start = readid*len;
		int start2 = RoundDown(&start);
		int shift_count = start - start2;
		mram_read(reads+(start2), read_cache, 160);
		// printf("%d read = %.150s\n", readid, read_cache+shift_count);
		
		// start alignment
		KmerIterator kit;
		KmerIteratornew(&kit, read_cache+shift_count);
		KmerIterator kit_end;
		KmerIteratornew_NULL(&kit_end);

		for (int i = 0;  !KmerIterator_cmp(&kit, &kit_end); ++i,KmerIterator_move(&kit)) {
			size_tt h = hash(&kit.kmer) & (size_-1);
			// toString(&kit.kmer);
			// printf("h %llu \n", h);
			size_tt find;
			for (;; h = (h+1!=size_ ? h+1 : 0)) {
				// printf("h %llu %llu\n", h, size_);
				__dma_aligned uint64_t table_kmer_cache;
				__mram_ptr void *target_addr = (table+(h*gap));
				//printf("addr %p\n", target_addr);
				mram_read(target_addr, &table_kmer_cache, sizeof(uint64_t));
				
				if (table_kmer_cache == EMPTY_KMER) {
					// empty slot, not in table
					// printf("empty %lu vs %lu\n", table_kmer_cache, EMPTY_KMER);
					find = size_;
					break;
				} else if (table_kmer_cache == kit.kmer.longs[0]) {
					// same key, found
					// printf("kmer %lu vs %lu\n", table_kmer_cache, kit.kmer.longs[0]);
					find = h;
					break;
				} 
			}
			if (find != size_) {
				__dma_aligned uint64_t table_T_cache[8];
				mram_read(table+(find*gap)+1, &table_T_cache, sizeof(uint64_t)*(gap-1));
				for(int chunk = 0; chunk < gap-1; chunk++) {
					// printf("%d chunk %lu\n", chunk, table_T_cache[chunk]);
					for(int t_count = 3; t_count >= 0; t_count--) {
						int16_t v = (table_T_cache[chunk] >> (16*t_count));
						// printf("%d t_count %d\n", t_count, v);
						if(v == -1) 
							break;
						uint8_t add = 1;
						for(int l = 0; l < t.len; l++) {
							if(t.T[l] == v) {
								add = 0;
								break;
							}
						}
						if(add)
							t.T[t.len++] = v;
					}
				}
			}
		}	
		mram_write(&t, &result[readid], sizeof(dpu_result));
		for(int i = 0; i < t.len; i++) {
			printf("%d read match %d\n", readid, t.T[i]);
		}
	}
	return 0;
}