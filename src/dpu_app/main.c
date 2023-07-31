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
__host size_tt size_;

// runing info (k, t_max, kmer_max)
__host int k;
__host int kmer_max;
__host int t_max;

//STDOUT_BUFFER_INIT(40960)

int RoundDown(int* a)
{
    return *(a) & (-8);
}

int RoundUp(int* a)
{
    return (*a + 3) & ~0x03;
}

void intersection(int16_t *a, int16_t a_len, dpu_result *t, uint8_t *first) {
	if(first) {
		memcpy(t->T, a, a_len*sizeof(int16_t));
		t->kmer = 1;
		assert(a_len <= T_LEN);
		t->len = a_len;
		first = 0;
	}
	else {
		int16_t new[T_LEN];
		int16_t  new_len = 0;
		int16_t i = 0, j = 0;
		while(i < t->len && j < a_len) {
			if(a[j] == t->T[i]) {
				new[new_len ++] = a[j];
				j ++;
				i ++;
			}
			else if(a[j] < t->T[i]) {
				j ++;
			}
			else {
				i ++;
			}
		}
		if(new_len != 0)
			t->kmer ++;
		memcpy(t->T, new, new_len*sizeof(int16_t));
		t->len = new_len;
		assert(new_len <= T_LEN);
	}
}

int main(){

	// init
	buddy_init(4096);
	unsigned int tasklet_id = me();
	int gap = t_max/4 + 1;

	for(int readid = tasklet_id; readid < PACKET_SIZE; readid+=NR_TASKLETS){
		// init result
		__dma_aligned dpu_result *t = buddy_alloc(sizeof(dpu_result));
		memset(t, 0, sizeof(dpu_result));
		// read out the read from mram
		int len = READ_LEN;
		__dma_aligned char read_cache[160];
		int start = readid*len;
		int start2 = RoundDown(&start);
		int shift_count = start - start2;
		mram_read(reads+(start2), read_cache, 160);
		// printf("%d read =\n", readid);
		
		// start alignment
		KmerIterator kit;
		KmerIteratornew(&kit, read_cache+shift_count);
		read_cache[shift_count+READ_LEN] = '\0';
		KmerIterator kit_end;
		KmerIteratornew_NULL(&kit_end);

		for (int i = 0;  !KmerIterator_cmp(&kit, &kit_end); ++i,KmerIterator_move(&kit)) {
			int16_t t_per_kmer[30];
			int16_t t_per_kmer_len = 0;
			size_tt h = hash(&kit.kmer) & (size_-1);
			// toString(&kit.kmer);
			// printf("h %llu \n", h);
			size_tt find;
			uint8_t first = 1;
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
						int l = 0;
						for(; l < t_per_kmer_len; l++) {
							if(t_per_kmer[l] == v) {
								add = 0;
								break;
							} 
							else if (t_per_kmer[l] > v)
								break;
						}
						if(add) {
							int16_t t_tmp[T_LEN];
							if(l != 0)
								memcpy(t_tmp, t_per_kmer, l*sizeof(int16_t));
							t_per_kmer[l] = v;
							if((t_per_kmer_len-l) != 0)
								memcpy(t_tmp+l+1, t_per_kmer+l, (t_per_kmer_len-l)*sizeof(int16_t));
							t_per_kmer_len++;
						}
					}
				}
				// intersection
				intersection(t_per_kmer, t_per_kmer_len, t, &first);
			}
		}	
		assert(t->len >= 0);
		assert(t->len <= T_LEN);
		// printf("%d read match %d kmer with %d len t: ", readid, t->kmer, t->len);
		for(int i = 0; i < t->len; i++) {
			assert(t->T[i] >= 0);
			//printf("%d ", t->T[i]);
		}
		// printf("\n");

		__mram_ptr void *target_addr = (result+readid);
		mram_write(t, target_addr, sizeof(dpu_result));
		
		buddy_free(t);
	}
	return 0;
}