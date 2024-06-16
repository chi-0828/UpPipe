#include "dpu_def.h"
#include "hash.h"
#include <mram.h>
#include <assert.h>
#include <alloc.h>
#include "defs.h"
#include "mutex.h"
#include <barrier.h>

// about 1MB read packet (from host)
__mram_noinit char reads[PACKET_CAPACITY];

// result (to host)
__mram_noinit dpu_result results[PACKET_SIZE];

// 60 MB hash table (from hsot)
__mram_noinit uint64_t table_key[MAX_TABLE_N];
__mram_noinit uint64_t table_value[MAX_TABLE_N];

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
void sort_by_count(dpu_result *t, int16_t *count_map) {
	int16_t new_len = 0;
	int16_t max = 0, max_i = 0;
	for(int16_t i = 0; i < COUNT_LEN; i ++) {
		if(count_map[i] > max) {
			max = count_map[i];
			max_i = i;
		}
	}
	if(max != 0) {
		for(int16_t i = 0; i < COUNT_LEN; i ++) {
			if(count_map[i] == max) {
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

	for(int readid = tasklet_id; readid < PACKET_SIZE; readid+=NR_TASKLETS){
		// init result for this read
		__dma_aligned dpu_result t;
		memset(&t, 0, sizeof(dpu_result));
		int16_t count_map[COUNT_LEN];
		memset(&count_map, 0, sizeof(int16_t)*COUNT_LEN);

		// get the RNA read from mram to WRAM
		int len = READ_LEN;
		__dma_aligned char read_cache[WRAM_READ_LEN];
		int start = readid*len;
		int start2 = RoundDown(&start);
		int shift_count = start - start2;
		mram_read(reads+(start2), read_cache, WRAM_READ_LEN);
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
			size_tt start = (h == 0 ? size_ - 1 : h -1);
			
			for (;; h = (h+1!=size_ ? h+1 : 0)) {
				assert(h != start);
				// printf("h %llu %llu\n", h, size_);
				int16_t wram_index = 0, flag = 0;
				// pre-fetch more data into WRAM for reducing data movement between WRAM and MRAM
				__dma_aligned uint64_t table_kmer_cache[WRAM_PREFETCH_SIZE];
				__mram_ptr void *target_addr = (table_key + h);
				mram_read(target_addr, &table_kmer_cache, sizeof(uint64_t)*WRAM_PREFETCH_SIZE);
				while (wram_index < WRAM_PREFETCH_SIZE) {
					// __mram_ptr void *target_addr = (table_key + h + wram_index);
					// mram_read(target_addr, &table_kmer_cache, sizeof(uint64_t));
					// printf("kmer %lu vs %lu\n", table_kmer_cache, kit.kmer.longs[0]);

					if (table_kmer_cache[wram_index] == EMPTY_KMER) {
						// empty slot, not in table
						flag = 1;
						break;
					} else if (table_kmer_cache[wram_index] == kit.kmer.longs[0]) {
						// same key, found
						// get information in hash table
						for(int chunk = 0; chunk < t_max/4; chunk++) {
							// reuse table_kmer_cache for reading table
							// need value of the found key
							mram_read(table_value+h+chunk, &table_kmer_cache, sizeof(uint64_t));
							
							int16_t t1 = (int16_t)((((uint64_t)table_kmer_cache[0])>>(0))&(0xFFFF));
							int16_t t2 = (int16_t)((((uint64_t)table_kmer_cache[0])>>(16))&(0xFFFF));
							int16_t t3 = (int16_t)((((uint64_t)table_kmer_cache[0])>>(32))&(0xFFFF));
							int16_t t4 = (int16_t)((((uint64_t)table_kmer_cache[0])>>(48))&(0xFFFF));
					
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
							assert(arr_id < COUNT_LEN);
							count_map[arr_id] ++;
						}
						
						flag = 1;
						break;
					}
					wram_index ++;
					h++;
				}
				if(flag)
					break;
			}
		}	
		
		sort_by_count(&t, count_map);
		// move the result to MRAM
		mram_write(&t, results+(readid), sizeof(dpu_result));

	}
	return 0;
}
