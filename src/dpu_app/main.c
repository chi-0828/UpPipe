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
// __mram_noinit int16_t t_arr[T_LEN*PACKET_SIZE];
// __host int32_t kmer_arr[PACKET_SIZE];
// __host int t_len_arr[PACKET_SIZE];

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

void sort_by_count(dpu_result *t) {
	int16_t new_len = 0;
	int16_t max = 0, max_i = 0;
	for(int16_t i = 0; i < T_LEN; i ++) {
		// if(me() == 1)
			// printf("%d ", t->T[i]);
		if(t->T[i] > max) {
			max = t->T[i];
			max_i = i;
		}
	}
	// if(me() == 1)
		// printf("\n");
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
	// buddy_init(512);
	unsigned int tasklet_id = me();
	// int read_chunk = PACKET_SIZE / NR_TASKLETS;
	// int read_chunk_remain = PACKET_SIZE % NR_TASKLETS;
	// int read_chunk_size = (tasklet_id < read_chunk_remain) ? read_chunk + 1 : read_chunk;
	// int tasklet_start = tasklet_id*
	int gap = t_max/4 + 1;

	for(int readid = tasklet_id; readid < PACKET_SIZE; readid+=NR_TASKLETS){
		// init result
		// __dma_aligned dpu_result *t = buddy_alloc(sizeof(dpu_result));
		// __dma_aligned dpu_result *t = mem_alloc(sizeof(dpu_result));
		__dma_aligned dpu_result t;
		// __dma_aligned dpu_result t;
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

		for (int i = 0;  !KmerIterator_cmp(&kit, &kit_end); ++i,KmerIterator_move(&kit)) {
			// int16_t *t_per_kmer = buddy_alloc(T_LEN*sizeof(int16_t));
			int16_t t_per_kmer[T_LEN];
			// int16_t *t_per_kmer = mem_alloc(T_LEN*sizeof(int16_t));
			int16_t t_per_kmer_len = 0;
			// size_tt h = hash(&kit.kmer) & (size_-1);
			size_tt h = kit.kmer.longs[0] & (size_-1);
			// toString(&kit.kmer);
			// printf("h %llu \n", h);
			size_tt find = size_, start = (h == 0 ? size_ - 1 : h -1);
			for (;; h = (h+1!=size_ ? h+1 : 0)) {
				assert(h != start);
				// printf("h %llu %llu\n", h, size_);
				__dma_aligned uint64_t table_kmer_cache;
				__mram_ptr void *target_addr = (table+(h*gap));
				//printf("addr %p\n", target_addr);
				mram_read(target_addr, &table_kmer_cache, sizeof(uint64_t));
				// printf("kmer %lu vs %lu\n", table_kmer_cache, kit.kmer.longs[0]);
				// break;
				if (table_kmer_cache == EMPTY_KMER) {
					// empty slot, not in table
					//printf("empty %lu vs %lu\n", table_kmer_cache, EMPTY_KMER);
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
				// __dma_aligned uint64_t table_T_cache[T_LEN/4];
				// mram_read(table+(find*gap)+1, table_T_cache, sizeof(uint64_t)*(gap-1));
				// int16_t t_per_kmer_len2 = t_per_kmer_len;
				for(int chunk = 0; chunk < gap-1; chunk++) {
					// printf("%d chunk %lu\n", chunk, table_T_cache[chunk]);
					__dma_aligned uint64_t table_T_cache;
					mram_read(table+(find*gap)+1+chunk, &table_T_cache, sizeof(uint64_t));
					int16_t t1 = (int16_t)((((uint64_t)table_T_cache)>>(0))&(0xFFFF));
					int16_t t2 = (int16_t)((((uint64_t)table_T_cache)>>(16))&(0xFFFF));
					int16_t t3 = (int16_t)((((uint64_t)table_T_cache)>>(32))&(0xFFFF));
					int16_t t4 = (int16_t)((((uint64_t)table_T_cache)>>(48))&(0xFFFF));
			
					// if( tasklet_id == 1 )
					// 	printf("%d %d %d %d\n", t4, t3, t2, t1);
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

				// intersection
				for(int16_t tid = 0; tid < t_per_kmer_len; tid ++) {
					int16_t arr_id = t_per_kmer[tid] - (int16_t)first_tid;
					t.T[arr_id] ++;
				}
				// buddy_free(t_per_kmer);
			}
		}	
		
		sort_by_count(&t);
		// if(me() == 1) {
		// 	printf("%d read match %d kmer with %d len t: ", readid, t.kmer, t.len);
		// 	for(int i = 0; i < t.len; i++) {
		// 		assert(t.T[i] >= 0);
		// 		printf("%d ", t.T[i]);
		// 	}
		// 	printf("\n");
		// }
		mram_write(&t, results+(readid), sizeof(dpu_result));
		
		// buddy_free(t);
		// mem_reset();
	}
	return 0;
}