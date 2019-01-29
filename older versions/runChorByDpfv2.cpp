#include <stdlib.h>
#include <stdio.h>
#include <err.h>
#include <unistd.h>
#include <immintrin.h>
#include <string.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <time.h>
#include <stdint.h>
#include <inttypes.h>
#include <iostream>
#include <bitset>
#include <initializer_list>

#include <cstdio>
#include <cstring>

#include "dpfv5.h"
//#include <fcntl.h>
//#include <sysexits.h>
//#include <sys/types.h>
//#include <sys/stat.h>


template <size_t T>
struct word { 
  static constexpr size_t number_of_chunk_for_a_word = T%32==0?T/32:T/32+1;
  __m256i chunk_for_a_word[number_of_chunk_for_a_word];
};

template <size_t B, size_t W>
struct record {
  static constexpr size_t number_of_struct_word_in_a_record = W;
  struct word <B> words_in_record[number_of_struct_word_in_a_record];
};

int main(int argc, char *argv[]) {
    int const number_of_rows = (1ULL << 10);//constant
    int number_of_bytes_per_row = atoi(argv[2]);//1500
    int const number_of_bytes_per_word = 50;//constant
    printf("Rows: %i, Bytes per row: %i, bytes per word: %i\n", number_of_rows, number_of_bytes_per_row, number_of_bytes_per_word);

    int const number_of_words_per_row = 30;//number_of_bytes_per_row % number_of_bytes_per_word == 0 ? number_of_bytes_per_row / number_of_bytes_per_word : number_of_bytes_per_row / number_of_bytes_per_word + 1;;
    int number_of_bytes_in_file = number_of_rows*number_of_bytes_per_row;
    int number_of_256_in_a_word = number_of_bytes_per_word%32==0?number_of_bytes_per_word/32:number_of_bytes_per_word/32+1;
    
    struct word <number_of_bytes_per_word> *chunks;
    int err = posix_memalign((void**)&chunks, 32, number_of_rows*number_of_words_per_row*sizeof(struct word<number_of_bytes_per_word>));
    if(err) perror("Error in memalign");
    int fd = open(argv[4], O_RDONLY);
    if (fd == -1) perror("Error opening file for reading");

    void *temp = NULL;
    char *data;
    temp = mmap(NULL, number_of_bytes_in_file, PROT_READ,MAP_SHARED, fd, 0);
    if(temp == MAP_FAILED) perror("Error in mmap");
    data = (char *) temp;
    for(int i=0;i<number_of_words_per_row*number_of_rows;++i) {
      memcpy(&chunks[i],data+i*number_of_bytes_per_word,number_of_bytes_per_word);
    }
    munmap(data,number_of_bytes_in_file);
    close(fd);

    struct record <number_of_bytes_per_word, number_of_words_per_row> *records;
    records = (struct record <number_of_bytes_per_word, number_of_words_per_row> *) chunks;
    
    //Chor's protocol by dpf

    //common to client and servers
    using namespace dpf;
    AES_KEY aeskey;
    AES_set_encrypt_key(_mm_set_epi64x(597349, 121379), &aeskey);

    typedef bool leaf_type;
    const size_t nitems = number_of_rows;//(1ULL << 12);
    dpf::dpf_key<nitems, leaf_type> dpfkey[2];
    __m128i * s, * s1;
    // printf("output_length: %d\n", (int)dpf_key<nitems,leaf_type>::output_length);
    if (posix_memalign((void**)&s, sizeof(__m128i), dpf_key<nitems,leaf_type>::output_length * sizeof(__m128i)) != 0)
      printf("alloc failed\n");
    if (posix_memalign((void**)&s1, sizeof(__m128i), dpf_key<nitems,leaf_type>::output_length * sizeof(__m128i)) != 0)
      printf("alloc failed\n");
    uint8_t * t = (uint8_t*)malloc(dpf_key<nitems,leaf_type>::output_length);
    uint8_t * t1 = (uint8_t*)malloc(dpf_key<nitems,leaf_type>::output_length);

    int query_length = number_of_rows % 64 == 0 ? number_of_rows/64:number_of_rows/64+1;

    //client query construction
    const size_t item = 501;//row (index) to fetch
    dpf::gen(aeskey, item, dpfkey, (leaf_type)1);

    //server 1 query expansion and response construction
    dpf::evalfull(aeskey, dpfkey[0], s, t);
    
    uint64_t query[query_length];
    
    for(int i=0;i<query_length/2;i++) {
      query[i*2]  = s[i][0];
      query[i*2+1] = s[i][1];
    }
    struct record <number_of_bytes_per_word, number_of_words_per_row> result;
    memset(&result, 0, sizeof(struct record <number_of_bytes_per_word, number_of_words_per_row>));

    int record_count = 0;
    for (size_t k = 0; k < query_length; ++k) {
      uint64_t bitset = query[k];
      // printf("%d\n", record_count);
      while (bitset != 0) {

        const int nextbit = __builtin_ctzll(bitset);//trailing zero, i.e. the row number where is a 1
        if(k*64+nextbit>=number_of_rows){break;}
        for(int p=0;p<number_of_words_per_row;p++) {
          for(int q=0;q<number_of_256_in_a_word;q++) {
            result.words_in_record[p].chunk_for_a_word[q] = _mm256_xor_si256(result.words_in_record[p].chunk_for_a_word[q],records[k*64+nextbit].words_in_record[p].chunk_for_a_word[q]);
          }
        }
        // std::cout << "bitset: " << std::bitset<64>(bitset) << "\n";
        // printf("nextbit: %d\n", nextbit );
        
        bitset ^= bitset & -bitset;//vanishes LSB 1
        record_count++;
      }
    }
    //server 2 query expansion and response construction
    dpf::evalfull(aeskey, dpfkey[1], s1, t1);
    uint64_t query1[query_length];
    
    for(int i=0;i<query_length/2;i++) {
      query1[i*2]  = s1[i][0];
      query1[i*2+1] = s1[i][1];
    }
    struct record <number_of_bytes_per_word, number_of_words_per_row> result1;
    memset(&result1, 0, sizeof(struct record <number_of_bytes_per_word, number_of_words_per_row>));

    int record_count1 = 0;
    for (size_t k = 0; k < query_length; ++k) {
      uint64_t bitset = query1[k];
      // printf("%d\n", record_count1);
      while (bitset != 0) {

        const int nextbit = __builtin_ctzll(bitset);//trailing zero, i.e. the row number where is a 1
        if(k*64+nextbit>=number_of_rows){break;}
        for(int p=0;p<number_of_words_per_row;p++) {
          for(int q=0;q<number_of_256_in_a_word;q++) {
            result1.words_in_record[p].chunk_for_a_word[q] = _mm256_xor_si256(result1.words_in_record[p].chunk_for_a_word[q],records[k*64+nextbit].words_in_record[p].chunk_for_a_word[q]);
          }
        }
        // std::cout << "bitset: " << std::bitset<64>(bitset) << "\n";
        // printf("nextbit: %d\n", nextbit );
        
        bitset ^= bitset & -bitset;//vanishes LSB 1
        record_count1++;
      }
    }

    //client record reconstruction
    struct record <number_of_bytes_per_word, number_of_words_per_row> fetched_record;
    memset(&fetched_record, 0, sizeof(struct record <number_of_bytes_per_word, number_of_words_per_row>));

    for(int p=0;p<number_of_words_per_row;p++) {
      for(int q=0;q<number_of_256_in_a_word;q++) {
        fetched_record.words_in_record[p].chunk_for_a_word[q] = _mm256_xor_si256(result.words_in_record[p].chunk_for_a_word[q], result1.words_in_record[p].chunk_for_a_word[q]);
      }
    }

    //Correctness
    struct record <number_of_bytes_per_word, number_of_words_per_row> correct_ness;
    memset(&correct_ness, 0, sizeof(struct record <number_of_bytes_per_word, number_of_words_per_row>));
    bool correct_algo = true;
    for(int p=0;p<number_of_words_per_row;p++) {
      for(int q=0;q<number_of_256_in_a_word;q++) {
        correct_ness.words_in_record[p].chunk_for_a_word[q] = _mm256_xor_si256(fetched_record.words_in_record[p].chunk_for_a_word[q], records[item].words_in_record[p].chunk_for_a_word[q]);
        if(correct_ness.words_in_record[p].chunk_for_a_word[q][0]!=0 || correct_ness.words_in_record[p].chunk_for_a_word[q][1]!=0 || correct_ness.words_in_record[p].chunk_for_a_word[q][2]!=0 || correct_ness.words_in_record[p].chunk_for_a_word[q][3]!=0){
          correct_algo = false;
          break;
        }        
      }
    }
    if(correct_algo) printf("Correctness successful!\n");
    else printf("Incorrect result!\n");
    // for(int i = 0; i < dpf_key<nitems,leaf_type>::output_length; ++i){
    //   printf("%llx\n",  s[i][0] ^ s1[i][0]);
    //   printf("%llx\n",  s[i][1] ^ s1[i][1]);
    // }
    
    free(s);
    free(t);

    free(s1);
    free(t1);


    
    free(chunks);
    return 0;
}

/*void splice(const _m128i * mask, uint8_t nmasks, _m1024i out)
{
	const _m256i shuffle = mm256_setr_epi64x(0x0000000000000000,
		0x0101010101010101, 0x0202020202020202, 0x0303030303030303);
	const _m256i bit_mask = mm256_set1_epi64x(0x7fbfdfeff7fbfdfe);
	const _m256i ff = mm256_set1_epi64x(0xffffffffffffffff);
	const _m256i shift[] = { mm256_set1_epi8(1), _mm256_set1_epi8(2),
		mm256_set1_epi8(4), mm256_set1_epi8(8), _mm256_set1_epi8(16),
		mm256_set1_epi8(32), mm256_set1_epi8(64), _mm256_set1_epi8(128) };

	for (uint8_t i = 0; i < nmasks; ++i)
	{
		uint32_t * maski = (uint32_t*)&mask[i];
		for (int j = 0; j < 4; j++)
		{
            out[j] = _mm256_or_si256(out[j], // accumulate into result
			_mm256_and_si256(shift[i], // zero all but i'th bit of each byte
				_mm256_cmpeq_epi8( // 0xff if byte is 0xff; else 0x00
					_mm256_or_si256( // set bits not possibly set by shuffle
						_mm256_shuffle_epi8( // shuffle 32 bits into 32 bytes
							_mm256_set1_epi32(maski[j]), shuffle),
						bit_mask),
					ff)
				)
			);
		}
	}
}*/