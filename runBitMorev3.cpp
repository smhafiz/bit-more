#include "bitmorev3.h"
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/time.h>

template <size_t nbytes_per_word>
struct word
{ 
  static constexpr size_t len256 = std::ceil(nbytes_per_word / static_cast<double>(sizeof(__m256i)));
  __m256i m256[len256];
  inline word<nbytes_per_word> & operator^=(const word<nbytes_per_word> & other)
  {
    for (size_t i = 0; i < len256; ++i)
    {
      m256[i] = _mm256_xor_si256(m256[i], other.m256[i]);
    }
    return *this;
  }
};

template <size_t nbytes_per_word, size_t nwords_per_row>
struct record
{
private:
  word<nbytes_per_word> words[nwords_per_row];
public:
  inline constexpr word<nbytes_per_word> & operator[](const size_t j)
  {
    return words[j];
  }
};

int main(int argc, char *argv[])
{
  AES_KEY aeskey;
  AES_set_encrypt_key(_mm_set_epi64x(597349, 121379), &aeskey);

  constexpr size_t lognitems = 12;
  constexpr size_t nitems = (1ULL << lognitems);
  constexpr size_t nbytes_per_row = 1024*30;
  constexpr size_t nkeys = 4;
  constexpr size_t nservers = std::exp2(nkeys);
  constexpr size_t nwords_per_row = nservers-1;
  constexpr size_t nbytes_per_word = std::ceil(nbytes_per_row / static_cast<double>(nwords_per_row));
  constexpr size_t alloc_size = nitems * sizeof(record<nbytes_per_word, nwords_per_row>);

  record<nbytes_per_word, nwords_per_row> * records;
  int err = posix_memalign((void**)&records, sizeof(__m256i), alloc_size);
  if(err) perror("Error in memalign for records");

  //read from file
  // int fd = open("testDB", O_RDONLY);
  // if (fd == -1) perror("Error opening file for reading");
  // char *data = (char *) mmap(NULL, alloc_size, PROT_READ,MAP_SHARED, fd, 0);
  // if(data == MAP_FAILED) perror("Error in mmap");
  // for(size_t i=0;i<nitems;++i) {
  //   memcpy(&records[i],data+i*sizeof(record<nbytes_per_word, nwords_per_row>),sizeof(record<nbytes_per_word, nwords_per_row>));
  // }
  // munmap(data,alloc_size);
  // close(fd);
  
  arc4random_buf(records, alloc_size);

  constexpr size_t item = 1;
  __m128i ** s = (__m128i**)malloc(nkeys * sizeof(__m128i *));
  uint8_t ** t = (uint8_t**)malloc(nkeys * sizeof(uint8_t *));
  dpf::dpf_key<nitems> dpfkey[nkeys][2];
  //size_t sindex = arc4random_uniform(1ULL << nkeys);
  dpf::dpf_key<nitems> dpf_query0[nkeys];
  dpf::dpf_key<nitems> dpf_query1[nkeys];
  for (size_t i = 0; i < nkeys; ++i)
  {
    err = posix_memalign((void**)&s[i], sizeof(__m256i), dpf::dpf_key<nitems>::output_length * sizeof(__m128i));
    if(err) perror("Error in memalign for s");
    t[i] = (uint8_t *)malloc(dpf::dpf_key<nitems>::output_length * sizeof(uint8_t));
    //dpf_query[i] = dpfkey[i][(sindex >> i) & 1U];
  }
  struct timeval t1, t2;
  double elapsed_time_keyGen;
  gettimeofday(&t1, NULL);
  size_t perm = bitmore::gen<nkeys>(aeskey, item, dpfkey);
  gettimeofday(&t2, NULL);
  elapsed_time_keyGen = (t2.tv_sec - t1.tv_sec) * 1000.0;
  elapsed_time_keyGen += (t2.tv_usec - t1.tv_usec) / 1000.0;
  printf("Key Gen for %d key pairs: %lf\n", nkeys, elapsed_time_keyGen);

  for (size_t i = 0; i < nservers; ++i)
  {
    printf("%zx, %zu\n", perm, perm^i);
  }

  for (size_t i = 0; i < nkeys; ++i)
  {
    dpf_query0[i] = dpfkey[i][((perm^0) >> i) & 1U];
    dpf_query1[i] = dpfkey[i][((perm^(nservers-1)) >> i) & 1U];
  }

  printf("\n%lu\n", bitmore::eval<nkeys, nservers>(aeskey, dpf_query1, item));

  //all servers
  word<nbytes_per_word> * results = (struct word<nbytes_per_word> *)malloc(nservers*sizeof(struct word<nbytes_per_word>));
  memset(results,0,nservers*sizeof(struct word<nbytes_per_word>));
  double elapsed_time_key_expand_and_response = 0.0;
  //size_t sindex = arc4random_uniform(1ULL << nkeys);
  for (size_t sindex = 0; sindex < nservers; ++sindex)
  {
    dpf::dpf_key<nitems> dpf_query_as[nkeys];
    size_t effective_perm = ((perm^sindex));
    for (size_t i = 0; i < nkeys; ++i)
    { size_t tempo = ( effective_perm >> i) & 1U;
      printf("%zu, ", tempo);
      dpf_query_as[i] = dpfkey[i][tempo];
    }
    printf("\n");
    gettimeofday(&t1, NULL);
    uint8_t * expanded_query;
    err = posix_memalign((void**)&expanded_query, sizeof(__m256i), nitems * sizeof(uint8_t));
    if(err) perror("Error in memalign for query");
    // dpf::evalfull_bitmore<nkeys, nitems, bool>(aeskey, dpf_query_as, s, t, expanded_query);
    bitmore::evalfull<nkeys,nservers>(aeskey, dpf_query_as, s, t, expanded_query);
    for (size_t i = 0; i < nitems; ++i)
    {
      if(expanded_query[i] < nwords_per_row) results[sindex] ^= records[i][expanded_query[i]];
    }
    free(expanded_query);
    gettimeofday(&t2, NULL);
    elapsed_time_key_expand_and_response += (t2.tv_sec - t1.tv_sec) * 1000.0;
    elapsed_time_key_expand_and_response += (t2.tv_usec - t1.tv_usec) / 1000.0;
  }
  printf("Key expansion and response generation per server: %lf\n", elapsed_time_key_expand_and_response/nservers);
  //record reconstruction
  for (size_t i = 0; i < nwords_per_row; ++i)
  {
    results[i]^=results[nservers-1];
  }
  //check correctness
  bool correct_ness = true;
  for (size_t i = 0; i < nwords_per_row; ++i)
  {
    results[i] ^= records[item-1][i];
    if(results[i].m256[0][0]!=0x00) {correct_ness = false; printf("Incorrect!\n");break;}
  }
  if (correct_ness) printf("Correct Protocol!\n");

  for (size_t i = 0; i < nkeys; ++i)
  {
    free(s[i]);
    free(t[i]);
  }
  free(s);
  free(t);

  free(records);
  free(results);
  return 0;
}