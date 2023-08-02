#include "dpu_def.h"
#include "hash.h"
#include <mram.h>
#include <assert.h>
#include <alloc.h>
#include "defs.h"
#include "mutex.h"
#include <barrier.h>

// about 1MB read packet (from hsot)
__mram_noinit char reads[PACKET_CAPACITY];

// result (to host)
__mram_noinit dpu_result results[PACKET_SIZE];

// 60 MB hash table (from hsot)
__mram_noinit uint64_t table[MAX_table_n];

// read info (from hsot)
__host size_tt size_;

// runing info (k, t_max)
__host int32_t k;
__host int32_t t_max;
__host int32_t first_tid;

// STDOUT_BUFFER_INIT(1024)

int RoundDown(int* a)
{
    return *(a) & (-8);
}

int RoundUp(int* a)
{
    return (*a + 3) & ~0x03;
}

// get the most similar transcript
void sort_by_count(dpu_result *t) {
	int16_t new_len = 0;
	int16_t max = 0, max_i = 0;
	for(int16_t i = 0; i < T_LEN; i ++) {
		if(t->T[i] > max) {
			max = t->T[i];
			max_i = i;
		}
	}
	if(max != 0) {
		for(int16_t i = 0; i < T_LEN; i ++) {
			if(t->T[i] == max) {
				t->T[new_len] = i + (int16_t)first_tid;;
				new_len++;
			}
		}
	}
	t->kmer = max;
	t->len = new_len;
}

int main(){

	// init
	unsigned int tasklet_id = me();
	int gap = t_max/4 + 1;

	for(int readid = tasklet_id; readid < PACKET_SIZE; readid+=NR_TASKLETS){
		// init result for this read
		__dma_aligned dpu_result t;
		memset(&t, 0, sizeof(dpu_result));

		// get the RNA read from mram to WRAM
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
		read_cache[shift_count+READ_LEN] = '\0';
		KmerIterator kit_end;
		KmerIteratornew_NULL(&kit_end);

		// iterate all kmers in this read
		for (int i = 0;  !KmerIterator_cmp(&kit, &kit_end); ++i,KmerIterator_move(&kit)) {
			int16_t t_per_kmer[T_LEN];
			int16_t t_per_kmer_len = 0;

			// find the kmer in hash table
			size_tt h = hash(&kit.kmer) & (size_-1);

			size_tt find = size_, start = (h == 0 ? size_ - 1 : h -1);
			for (;; h = (h+1!=size_ ? h+1 : 0)) {
				assert(h != start);
				// printf("h %llu %llu\n", h, size_);
				__dma_aligned uint64_t table_kmer_cache;
				__mram_ptr void *target_addr = (table+(h*gap));
				mram_read(target_addr, &table_kmer_cache, sizeof(uint64_t));
				// printf("kmer %lu vs %lu\n", table_kmer_cache, kit.kmer.longs[0]);

				if (table_kmer_cache == EMPTY_KMER) {
					// empty slot, not in table
					find = size_;
					break;
				} else if (table_kmer_cache == kit.kmer.longs[0]) {
					// same key, found
					find = h;
					break;
				} 
			}
			// get information in hash table
			if (find != size_) {
				// TODO: read more in one time
				for(int chunk = 0; chunk < gap-1; chunk++) {
					__dma_aligned uint64_t table_T_cache;
					mram_read(table+(find*gap)+1+chunk, &table_T_cache, sizeof(uint64_t));
					int16_t t1 = (int16_t)((((uint64_t)table_T_cache)>>(0))&(0xFFFF));
					int16_t t2 = (int16_t)((((uint64_t)table_T_cache)>>(16))&(0xFFFF));
					int16_t t3 = (int16_t)((((uint64_t)table_T_cache)>>(32))&(0xFFFF));
					int16_t t4 = (int16_t)((((uint64_t)table_T_cache)>>(48))&(0xFFFF));
			
					if(t4 != -1) {
						t_per_kmer[t_per_kmer_len++] = t4;
						if(t3 != -1){
							t_per_kmer[t_per_kmer_len++] = t3;
							if(t2 != -1){
								t_per_kmer[t_per_kmer_len++] = t2;
								if(t1 != -1){
									t_per_kmer[t_per_kmer_len++] = t1;
								}
							}
						}
					}
					assert(t_per_kmer_len <= T_LEN);
				}

				// count the transcript id
				for(int16_t tid = 0; tid < t_per_kmer_len; tid ++) {
					int16_t arr_id = t_per_kmer[tid] - (int16_t)first_tid;
					t.T[arr_id] ++;
				}
			}
		}	
		
		sort_by_count(&t);
		// move the result to MRAM
		mram_write(&t, results+(readid), sizeof(dpu_result));

	}
	return 0;
}