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
// __mram_noinit int64_t result_2_read[MAX_PACKET_SIZE];
// __mram_noinit int64_t result_id[MAX_PACKET_SIZE];
// __mram_noinit int64_t result_pos[MAX_PACKET_SIZE];
__host int result_len[PACKET_SIZE];

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

// multi-tasklets
BARRIER_INIT(my_barrier, NR_TASKLETS);
int64_t result_id_tasklet[NR_TASKLETS][16];
int64_t result_pos_tasklet[NR_TASKLETS][16];
int32_t v_16len[NR_TASKLETS];
int32_t head;
int32_t matched_size_tasklet[NR_TASKLETS];

//STDOUT_BUFFER_INIT(40960)

int RoundDown(int* a)
{
    return *(a) & (-8);
}

int main(){

	// init
	unsigned int tasklet_id = me();

	// read out the hash table
	KmerHashTable kmap;
	set_empty(&(kmap.empty));
	kmap.size_ = size_;
	kmap.table_ptr = table;

	// rseult cahce
	//   int64_t *vint = result_id_tasklet[tasklet_id];
	//   int64_t *vpos= result_pos_tasklet[tasklet_id];
	//   v_16len[tasklet_id] = 0;
	//   int v_matched = 0;
	//   int v_len = v_16len[tasklet_id];

	for(int readid = tasklet_id; readid < PACKET_SIZE; readid+=NR_TASKLETS){

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

			size_tt h = hash(&kit.kmer) & (kmap.size_-1);
			//printf("h %llu \n", h);
			//h += round_hash*MAX_table_n;
			size_tt find;
			for (;; h = (h+1!=kmap.size_ ? h+1 : 0)) {
				if(h+1!=kmap.size_ -1)
					break;
				__dma_aligned uint64_t table_kmer_cache;
				
				int gap = t_max/4 + 1;
				__mram_ptr void *target_addr = (kmap.table_ptr+(h*gap));
				//printf("addr %p\n", target_addr);p 
				mram_read(target_addr, &table_kmer_cache, sizeof(uint64_t));
				printf("kmer %lu vs %lu\n", table_kmer_cache, kmap.empty.longs[0]);
				
				if (table_kmer_cache == kmap.empty.longs[0]) {
					//*matched_kmer = table_kmer_cache[0];
					// empty slot, not in table
					find = kmap.size_;
				} else if (table_kmer_cache == kit.kmer.longs[0]) {
					//*matched_kmer = table_kmer_cache[0];
					// same key, found
					find = h;
				} // if it is deleted, we still have to continue
			}
			// if (finded_key != kmap.size_) {
      
			// 	// // read table to wram
			// 	// __dma_aligned int64_t table_int_cache[1];
			// 	// __mram_ptr void *target_addr = kmap.table_int_ptr+finded_key;
			// 	// mram_read(target_addr, table_int_cache, sizeof(int64_t));
			// 	// //printf("=> contig %lld ", table_int_cache[0]);

			// 	// // prevent error or bug
			// 	// assert(table_int_cache[0] != -1);
			// }
		}
		// // iterate each read
		// v_matched = match(readid, read_cache, shift_count, len, &kmap, vint, vpos, tasklet_id, result_id, result_pos);

		// // get matched count of this read
		// matched_count_read = (matched_size_tasklet[tasklet_id] - matched_count_read);
		// result_len[readid] = matched_count_read;
	}
	return 0;
}