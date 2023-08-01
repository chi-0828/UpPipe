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
__mram_noinit dpu_result result[PACKET_SIZE];
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
	int16_t new[T_LEN];
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
				new[new_len] = i + (int16_t)first_tid;;
				new_len++;
			}
		}
	}
	t->kmer = max;
	t->len = new_len;
	memcpy(t->T, new, sizeof(int16_t)*new_len);
}

void sort_by_id(int16_t arr[], int16_t arr_count[], int16_t n) {
	int16_t i, key, j;
	int16_t key2;
	// sort
    for (i = 1; i < n; i++) {
        key = arr[i];
		key2 = arr_count[i];
        j = i - 1;
        while (j >= 0 && arr[j] > key) {
            arr[j + 1] = arr[j];
			arr_count[j + 1] = arr_count[j];
            j = j - 1;
        }
        arr[j + 1] = key;
		arr_count[j + 1] = key2;
    }
}

void intersection(int16_t a[], int16_t a_len, dpu_result *t) {

	for(int16_t i = 0; i < a_len; i ++) {
		int16_t arr_id = a[i] - (int16_t)first_tid;
		t->T[arr_id] ++;
	}

}

int main(){

	// init
	buddy_init(4096);
	unsigned int tasklet_id = me();
	int gap = t_max/4 + 1;

	for(int readid = tasklet_id; readid < PACKET_SIZE; readid+=NR_TASKLETS){
		// init result
		//__dma_aligned dpu_result *t = buddy_alloc(sizeof(dpu_result));
		__dma_aligned dpu_result t;
		memset(&t, 0, sizeof(dpu_result));
		// get the RNA read from mram to WRAM
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
			int16_t *t_per_kmer = buddy_alloc(T_LEN*sizeof(int16_t));
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
			
					// printf("%d t_count %d\n", t_count, v);
					if(t4 != -1) {
						t_per_kmer[t_per_kmer_len++] = t4;
					}
					if(t3 != -1){
						t_per_kmer[t_per_kmer_len++] = t3;
					}
					if(t2 != -1){
						t_per_kmer[t_per_kmer_len++] = t2;
					}
					if(t1 != -1){
						t_per_kmer[t_per_kmer_len++] = t1;
					}
					assert(t_per_kmer_len <= T_LEN);
				}
				// intersection
				intersection(t_per_kmer, t_per_kmer_len, &t);
				buddy_free(t_per_kmer);
			}
		}	
		
		sort_by_count(&t);

		// printf("%d read match %d kmer with %d len t: ", readid, t.kmer, t.len);
		// for(int i = 0; i < t.len; i++) {
		// 	assert(t.T[i] >= 0);
		// 	printf("%d ", t.T[i]);
		// }
		// printf("\n");

		__mram_ptr void *target_addr = (result+readid);
		mram_write(&t, target_addr, sizeof(dpu_result));
		
		// buddy_free(t);
	}
	return 0;
}