#include "dpf.h"

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

uint8_t reverse_a_number(uint8_t x, size_t limit){
  uint8_t y = 0x00;
  for (int i = 0; i < limit; ++i)
  {
    y += (x>>i& 1U)*std::exp2(limit-1-i);
  }
  //printf("%x->%x\n", x, y);
  return y;
}

int main(int argc, char *argv[])
{
  AES_KEY aeskey;
  AES_set_encrypt_key(_mm_set_epi64x(597349, 121379), &aeskey);

  constexpr size_t nitems = (1ULL << 10);//number of rows
  constexpr size_t nbytes_per_row = 1024*15;//b bits
  constexpr size_t nkeys = 6;//L bits
  constexpr size_t nservers = 64;//ell
  constexpr size_t nwords_per_row = nservers-1;//s = ell - 1
  constexpr size_t nbytes_per_word = std::ceil(nbytes_per_row / static_cast<double>(nwords_per_row));

  //printf("Rows: %lu, Bytes per row: %lu, bytes per word: %lu\n", nitems, nbytes_per_row, nbytes_per_word);

  constexpr size_t alloc_size = nitems * sizeof(record<nbytes_per_word, nwords_per_row>);//DB size; r * s * w 

  record<nbytes_per_word, nwords_per_row> * records;
  posix_memalign((void**)&records, sizeof(__m256i), alloc_size);
  arc4random_buf(records, alloc_size);

  constexpr size_t item = 501;
  __m128i ** s = (__m128i**)malloc(nkeys * sizeof(__m128i *));
  uint8_t ** t = (uint8_t**)malloc(nkeys * sizeof(uint8_t *));

  dpf::dpf_key<nitems> dpfkey[nkeys][2];

  //size_t sindex = arc4random_uniform(1ULL << nkeys);

  dpf::dpf_key<nitems> dpf_query[nkeys];
  for (size_t i = 0; i < nkeys; ++i)
  {
    posix_memalign((void**)&s[i], sizeof(__m256i), dpf::dpf_key<nitems>::output_length * sizeof(__m128i));
    t[i] = (uint8_t *)malloc(dpf::dpf_key<nitems>::output_length * sizeof(uint8_t));
    dpf::gen(aeskey, item, dpfkey[i]);
    if (eval(aeskey,dpfkey[i][0],item))
    {
      dpf::dpf_key<nitems> temp = dpfkey[i][0];
      dpfkey[i][0] = dpfkey[i][1];
      dpfkey[i][1] = temp;
    }
    //dpf_query[i] = dpfkey[i][(sindex >> i) & 1U];
  }
  //all servers response generation
  word<nbytes_per_word> * result = (struct word<nbytes_per_word> *)malloc(nservers*sizeof(struct word<nbytes_per_word>));
  memset(result,0,nservers*sizeof(struct word<nbytes_per_word>));
  //size_t sindex = arc4random_uniform(1ULL << nkeys);
  for (size_t sindex = 0; sindex < nservers; ++sindex)
  {
    dpf::dpf_key<nitems> dpf_query[nkeys];
    for (size_t i = 0; i < nkeys; ++i)
    {
      dpf_query[i] = dpfkey[i][(sindex >> i) & 1U];
    }

    uint8_t * expanded_query;
    posix_memalign((void**)&expanded_query, sizeof(__m256i), nitems * sizeof(uint8_t));
    dpf::evalfull_bitmore<nkeys, nitems, bool>(aeskey, dpf_query, s, t, expanded_query);
    
    for (size_t i = 0; i < nitems; ++i)
    {
      if(expanded_query[i] < nwords_per_row) result[sindex] ^= records[i][expanded_query[i]];
    }

    free(expanded_query);

  }

  // uint8_t keymap = 0;
  // for (size_t i = 0; i < nkeys; ++i)
  // {
  //   //printf("%u\n", eval(aeskey,dpfkey[i][0],item));
  //   keymap += (std::exp2(i)) * eval(aeskey,dpfkey[i][0],item);
  // }
  //printf("%u\n", keymap);
  //record reconstruction
  for (size_t i = 0; i < nwords_per_row; ++i)
  {
    result[i]^=result[nservers-1];
  }
  //check correctness
  for (size_t i = 0; i < nwords_per_row; ++i)
  {
    uint8_t tmp = reverse_a_number(i,nkeys);
    result[tmp] ^= records[item][i];
    printf("%x->%x\n", i, tmp);
    printf("%llu\n",result[tmp].m256[0][0]);

  }
  for (size_t i = 0; i < nkeys; ++i)
  {
    free(s[i]);
    free(t[i]);
  }
  free(s);
  free(t);

  free(records);
  free(result);

  return 0;
}