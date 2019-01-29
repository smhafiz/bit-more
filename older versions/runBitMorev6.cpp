#include "bitmorev5.h"
#include <algorithm>

template <int begin, int end, typename lambda>
inline void static_for(lambda const & f)
{
  if constexpr (begin < end)
  {
    f(std::integral_constant<int, begin>{});
    static_for<begin + 1, end>(f);
  }
}

double getTimeElapsed(struct timeval end, struct timeval start)
{
  return ((end.tv_sec - start.tv_sec)*1000.00) + ((end.tv_usec - start.tv_usec) / 1000.00);
  //return ((end.tv_sec - start.tv_sec)*1000000.00) + ((end.tv_usec - start.tv_usec));
}

int main(int argc, char * argv[])
{
  int ntrials = 3, nmetrics = 4, number_of_trials_pruned=0;
  struct timeval tvalBefore, tvalAfter;
  AES_KEY aeskey;
  AES_set_encrypt_key(_mm_set_epi64x(597349, 121379), &aeskey);

  constexpr size_t nitems = (1ULL << 12);
  constexpr size_t nbytes_per_row = 1024*127;

  constexpr size_t item = 101;

  void * database_;
  posix_memalign((void**)&database_, sizeof(__m256i), nitems * nbytes_per_row * 1.5);
  arc4random_buf(database_, nitems * nbytes_per_row * 1.5);

  uint8_t * expanded_query;
  posix_memalign((void**)&expanded_query, sizeof(__m256i), nitems * sizeof(uint8_t));

  constexpr std::array<uint8_t, 7> Nservers = {2, 4, 8, 16, 32, 64, 128};//{ 3, 6, 10, 27, 55, 112, 127 };
  auto benchmark = [&](auto i)
    {
      constexpr uint8_t nservers = std::get<i.value>(Nservers);
      constexpr size_t nwords_per_row = nservers-1;
      constexpr size_t nbytes_per_word = std::ceil(nbytes_per_row / static_cast<double>(nwords_per_row));
      constexpr size_t soundness = 0;
      constexpr size_t nkeys = std::log2(nservers)+soundness;

      for(int trials=0;trials<ntrials;++trials)
      {
        __m128i ** s = (__m128i**)malloc(nkeys * sizeof(__m128i *));
        uint8_t ** t = (uint8_t**)malloc(nkeys * sizeof(uint8_t *));
        for (size_t i = 0; i < nkeys; ++i)
        {
          posix_memalign((void**)&s[i], sizeof(__m256i), dpf::dpf_key<nitems>::output_length * sizeof(__m128i));
          t[i] = (uint8_t *)malloc(dpf::dpf_key<nitems>::output_length * sizeof(uint8_t));
        }

        typedef bitmore::record<nbytes_per_word, nwords_per_row> record;
        typedef bitmore::query<nitems, nkeys, nservers> query;

        record * database = reinterpret_cast<record *>(database_);
        //Client: query generation
        query q(aeskey, item);

        //Query expansion at s+1 th server
        auto Q = q.get_request(nwords_per_row);
        bitmore::evalfull<nkeys,nitems,nservers>(aeskey, Q, s, t, expanded_query);

        //response word from s+1 th server
        bitmore::word<nbytes_per_word> qtemplate = { 0 };
        for (size_t i = 0; i < nitems; ++i)
        {
          if (expanded_query[i] != nwords_per_row) qtemplate ^= database[i][expanded_query[i]];
        }

        //Query expansion and response word generation at the rest s number of servers, 0...s-1
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

        for (size_t i = 0; i < nwords_per_row; ++i)
        {
          result[i] ^= database[item][i];
        }

        bool correct_ness = true;
        for (size_t i = 0; i < nwords_per_row; ++i)
        {
          if(result[i].m256[0][0]!=0x00) {correct_ness = false; printf("%zu->Incorrect!\n",nservers);break;}
        }
        if (correct_ness) printf("%zu->Correct protocol!\n",nservers);

        for (size_t i = 0; i < nkeys; ++i)
        {
          free(s[i]);
          free(t[i]);
        }
        free(s);
        free(t);

      }
          };
  static_for<0, 7>(benchmark);


  free(database_);
  free(expanded_query);

  return 0;
}