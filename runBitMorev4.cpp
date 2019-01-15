#include "bitmorev4.h"

int main(int argc, char *argv[])
{
  AES_KEY aeskey;
  AES_set_encrypt_key(_mm_set_epi64x(597349, 121379), &aeskey);

  constexpr size_t nitems = (1ULL << 12);
  constexpr size_t nbytes_per_row = 1024*15;
  constexpr size_t nkeys = 1000;
  constexpr size_t nservers = 128;
  constexpr size_t nwords_per_row = nservers-1;
  constexpr size_t nbytes_per_word = std::ceil(nbytes_per_row / static_cast<double>(nwords_per_row));
  constexpr size_t item = 1;

  typedef dpf::dpf_key<nitems> dpf_key;
  typedef bitmore::record<nbytes_per_word, nwords_per_row> record;
  typedef bitmore::query<nitems, nkeys, nservers> query;

  record * database;
  posix_memalign((void**)&database, sizeof(__m256i), nitems * sizeof(record));
  arc4random_buf(database, nitems * sizeof(record));

  __m128i ** s = (__m128i**)malloc(nkeys * sizeof(__m128i *));
  uint8_t ** t = (uint8_t**)malloc(nkeys * sizeof(uint8_t *));
  for (size_t i = 0; i < nkeys; ++i)
  {
    posix_memalign((void**)&s[i], sizeof(__m256i), dpf_key::output_length * sizeof(__m128i));
    t[i] = (uint8_t *)malloc(dpf_key::output_length * sizeof(uint8_t));
  }
  query q(aeskey, item);

  uint8_t * expanded_query;
  posix_memalign((void**)&expanded_query, sizeof(__m256i), nitems * sizeof(uint8_t));

  auto Q = q.get_request(nwords_per_row);
  bitmore::evalfull<nkeys,nitems,nservers>(aeskey, Q, s, t, expanded_query);

  bitmore::word<nbytes_per_word> qtemplate = { 0 };
  for (size_t i = 0; i < nitems; ++i)
  {
    if (expanded_query[i] != nwords_per_row) qtemplate ^= database[i][expanded_query[i]];
  }

  record result = { qtemplate };
  for (size_t k = 0; k < nservers-1; ++k)
  {
    Q = q.get_request(k);
    bitmore::evalfull<nkeys,nitems,nservers>(aeskey, Q, s, t, expanded_query);
    for (size_t i = 0; i < nitems; ++i)
    {
      if (expanded_query[i] != nwords_per_row) result[k] ^= database[i][expanded_query[i]];
    }
  }

  printf("%llu=?=%llu\n", result[0].m256[0][0], database[item][0].m256[0][0]);
  printf("%llu=?=%llu\n", result[1].m256[0][0], database[item][1].m256[0][0]);
  printf("%llu=?=%llu\n", result[2].m256[0][0], database[item][2].m256[0][0]);

  free(database);
  free(expanded_query);

  return 0;
}