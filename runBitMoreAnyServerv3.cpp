#include "bitmorev3.h"

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

  constexpr size_t nitems = (1ULL << 12);
  constexpr size_t nbytes_per_row = 1024*15;
  
  constexpr size_t nservers = 7;
  constexpr size_t soundness = 8;
  constexpr size_t nkeys = std::ceil(std::log2(nservers))+soundness;
  constexpr size_t nwords_per_row = nservers-1;
  constexpr size_t nbytes_per_word = std::ceil(nbytes_per_row / static_cast<double>(nwords_per_row));

  //printf("Rows: %lu, Bytes per row: %lu, bytes per word: %lu\n", nitems, nbytes_per_row, nbytes_per_word);

  constexpr size_t alloc_size = nitems * sizeof(record<nbytes_per_word, nwords_per_row>);

  record<nbytes_per_word, nwords_per_row> * records;
  posix_memalign((void**)&records, sizeof(__m256i), alloc_size);
  arc4random_buf(records, alloc_size);

  constexpr size_t item = 1;
  __m128i ** s = (__m128i**)malloc(nkeys * sizeof(__m128i *));
  uint8_t ** t = (uint8_t**)malloc(nkeys * sizeof(uint8_t *));
  dpf::dpf_key<nitems> dpfkey[nkeys][2];
  //size_t sindex = arc4random_uniform(1ULL << nkeys);
  // dpf::dpf_key<nitems> dpf_query0[nkeys];
  // dpf::dpf_key<nitems> dpf_query1[nkeys];
  for (size_t i = 0; i < nkeys; ++i)
  {
    posix_memalign((void**)&s[i], sizeof(__m256i), dpf::dpf_key<nitems>::output_length * sizeof(__m128i));
    t[i] = (uint8_t *)malloc(dpf::dpf_key<nitems>::output_length * sizeof(uint8_t));
    //dpf_query[i] = dpfkey[i][(sindex >> i) & 1U];
  }
  size_t perm = bitmore::gen<nkeys>(aeskey, item, dpfkey);
  //printf("%u\n", sizeof(size_t));
  //printf("%zx\n", perm);
  for (size_t i = 0; i < nservers; ++i)
  {
    printf("%zx, %zu\n", perm^i, (perm^i)%nservers);
  }

   //all servers
  word<nbytes_per_word> * results = (struct word<nbytes_per_word> *)malloc(nservers*sizeof(struct word<nbytes_per_word>));
  memset(results,0,nservers*sizeof(struct word<nbytes_per_word>));
  //size_t sindex = arc4random_uniform(1ULL << nkeys);
  for (size_t sindex = 0; sindex < nservers; ++sindex)
  {
    dpf::dpf_key<nitems> dpf_query_as[nkeys];
    size_t effective_perm = ((perm^sindex)%nservers);
    for (size_t i = 0; i < nkeys; ++i)
    { size_t tempo = ( effective_perm >> i) & 1ULL;
      printf("%u, ", tempo);
      dpf_query_as[i] = dpfkey[i][tempo];
    }
    printf("\n\n");
    uint8_t * expanded_query;
    posix_memalign((void**)&expanded_query, sizeof(__m256i), nitems * sizeof(uint8_t));
    // dpf::evalfull_bitmore<nkeys, nitems, bool>(aeskey, dpf_query_as, s, t, expanded_query);
    bitmore::evalfull<nkeys,nservers>(aeskey, dpf_query_as, s, t, expanded_query);
    for (size_t i = 0; i < nitems; ++i)
    {
      if(expanded_query[i] < nwords_per_row) results[sindex] ^= records[i][expanded_query[i]];
    }

    free(expanded_query);

  }

  //record reconstruction
  for (size_t i = 0; i < nwords_per_row; ++i)
  {
    results[i]^=results[nservers-1];
  }
  //check correctness
  bool correct_ness = true;
  for (size_t i = 0; i < nwords_per_row; ++i)
  {
    results[i] ^= records[item][i];
    //printf("Word[%u]->%llu \n",i,results[i].m256[0][0]);
    if(results[i].m256[0][0]!=0x00) {correct_ness = false; printf("Incorrect!\n");break;}
  }
  if (correct_ness) printf("Correct protocol!\n");

  for (size_t i = 0; i < nkeys; ++i)
  {
    free(s[i]);
    free(t[i]);
  }
  free(s);
  free(t);

  free(records);
  // free(expanded_query0);
  // free(expanded_query1);
  free(results);
  return 0;
}