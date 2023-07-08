#include "dpu_def.h"
#include "hash.h"
#include <mram.h>
#include <assert.h>
#include "defs.h"
#include "mutex.h"
#include <barrier.h>

// about 1MB read packet (from hsot)
__mram_noinit char reads[MAX_PACKET_SIZE];

// result (to host)
__mram_noinit int64_t result_2_read[MAX_PACKET_SIZE];
__mram_noinit int64_t result_id[MAX_PACKET_SIZE];
__mram_noinit int64_t result_pos[MAX_PACKET_SIZE];
__host int result_len[PACKET_SIZE];

// 60 MB hash table (from hsot)
__mram_noinit int64_t table_int[MAX_table_n];
__mram_noinit uint64_t table_kmer[MAX_table_n];

// read info (from hsot)
__host int32_t reads_len[PACKET_SIZE];
__host int32_t read_n;
__host size_tt size_;

// result info (to hsot)
__host int matched_size;

// runing info (k and what should dpu do)
__host dpu_args arg;

// for large k , we have to record result array head position
__host int32_t round_hash;
__host int32_t last_round;
__mram int64_t result_head;

// multi-tasklets
BARRIER_INIT(my_barrier, NR_TASKLETS);
MUTEX_INIT(my_mutex);
int64_t result_id_tasklet[NR_TASKLETS][16];
int64_t result_pos_tasklet[NR_TASKLETS][16];
int32_t v_16len[NR_TASKLETS];
int32_t head;
int32_t matched_size_tasklet[NR_TASKLETS];

//STDOUT_BUFFER_INIT(40960)
    
int match(int rid, const char *s, int start, int len, KmerHashTable* kmap, int64_t* v_int, int64_t* v_pos, unsigned int tid, __mram_ptr int64_t* result_id, __mram_ptr int64_t* result_pos) {

  int v_acc = 0;

  // kmer iterator
  KmerIterator kit;
  KmerIteratornew(&kit, s+start);
  KmerIterator kit_end;
  KmerIteratornew_NULL(&kit_end);

  // iterating reads
  for (int i = 0;  !KmerIterator_cmp(&kit, &kit_end); ++i,KmerIterator_move(&kit)) {
    // find the key
    Kmer key = rep(kit.kmer);
    //toString(&(kit.kmer));
    size_tt finded_key = find(&key, kmap);
    int pos = kit.p_;

    if (finded_key != kmap->size_) {
      
      // read table to wram
      __dma_aligned int64_t table_int_cache[1];
      __mram_ptr void *target_addr = kmap->table_int_ptr+finded_key;
      mram_read(target_addr, table_int_cache, sizeof(int64_t));
      //printf("=> contig %lld ", table_int_cache[0]);

      // prevent error or bug
      assert(table_int_cache[0] != -1);

      // push back to array
      
      v_int[v_16len[tid]] = table_int_cache[0];
      v_pos[v_16len[tid]] = pos;
      v_16len[tid]++;
      matched_size_tasklet[tid]++;
      if(v_16len[tid] == 16){
        // WRAM cache is full, store to MRAM
        // lock
        mutex_lock(my_mutex);
        // __dma_aligned int64_t *vint_tmp =  result_id_tasklet[tid];
        // __dma_aligned int64_t *vpos_tmp =  result_pos_tasklet[tid];
        mram_write(v_int, result_id + head, 128);
        mram_write(v_pos, result_pos + head, 128);
        __dma_aligned int64_t r_2_r_tmp[16];
        for(int ii = 0; ii < 16; ii++){
          r_2_r_tmp[ii] = (int64_t)rid;
        }
        // for(int ii = 0; ii < 16; ii++){
        //   printf("%lld ", r_2_r_tmp[ii]);
        // }
        mram_write(r_2_r_tmp, result_2_read + head, 128);
        head += (16);
        matched_size = head;
        // unlock
        mutex_unlock(my_mutex);
        // reset cache
        v_16len[tid] = 0;
      }
      
    }
    else{
      
      //printf("not find ");
    }
    //printf("i(%d) %d\n", i, (len - k + 1));
    if(i == (len - arg.k + 1) -1)
      break;
  }

  // store the result from WRAM cache to MRAM
  // lock
  //mutex_unlock(my_mutex);
  if(v_16len[tid] > 0){
    mutex_lock(my_mutex);
    // __dma_aligned int64_t *vint_tmp =  result_id_tasklet[tid];
    // __dma_aligned int64_t *vpos_tmp =  result_pos_tasklet[tid];
    mram_write(v_int, result_id + head, v_16len[tid]*8);
    mram_write(v_pos, result_pos + head, v_16len[tid]*8);
    __dma_aligned int64_t r_2_r_tmp[16];
    for(int ii = 0; ii < 16; ii++){
      r_2_r_tmp[ii] = (int64_t)rid;
    }
    // for(int ii = 0; ii < 16; ii++){
    //   printf("%lld ", r_2_r_tmp[ii]);
    // }
    mram_write(r_2_r_tmp, result_2_read + head, v_16len[tid]*8);
    head += (v_16len[tid]);
    matched_size = head;
    // unlock
    mutex_unlock(my_mutex);
    // reset cache
    v_16len[tid] = 0;
  }

  //mutex_unlock(my_mutex);

  //printf("\n");
  return 1;
}

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
  kmap.table_int_ptr = table_int;
  kmap.table_kmer_ptr = table_kmer;

  // rseult cahce
//   int64_t *vint = result_id_tasklet[tasklet_id];
//   int64_t *vpos= result_pos_tasklet[tasklet_id];
//   v_16len[tasklet_id] = 0;
//   int v_matched = 0;
//   int v_len = v_16len[tasklet_id];

  for(int readid = tasklet_id ; readid < read_n; readid+=NR_TASKLETS){

    // read out the read from mram
    int len = reads_len[readid];
    __dma_aligned char read_cache[160];
    int start = readid*len;
    int start2 = RoundDown(&start);
    int shift_count = start - start2;
    mram_read(reads+(start), read_cache, 160);
	printf("%d read = %s\n", readid, read_cache);
    // // record the previous matched count 
    // int matched_count_read = matched_size_tasklet[tasklet_id];

    // // iterate each read
    // v_matched = match(readid, read_cache, shift_count, len, &kmap, vint, vpos, tasklet_id, result_id, result_pos);

    // // get matched count of this read
    // matched_count_read = (matched_size_tasklet[tasklet_id] - matched_count_read);
    // result_len[readid] = matched_count_read;
  }
  return 0;
}