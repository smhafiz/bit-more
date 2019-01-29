#ifndef BITMORE_H__
#define BITMORE_H__

#include "dpfv3.h"

namespace bitmore
{

using namespace dpf;

/*which server gives me zero, solve it
  find that perm
  nkeys bit string
  such that when it will xor with perm, reduce it mod nservers give zero
  that is the server for zeroth word
  xor perm with sever index and reduce it to mod nservers
*/
static const __m256i lo4 = _mm256_set1_epi8(0x0f);

template <size_t nkeys, size_t nitems>
size_t gen(AES_KEY & aeskey, size_t point, dpf_key<nitems> dpfkey[nkeys][2])
{
  size_t perm;
  memset(&perm,0,sizeof(size_t));
  for (size_t i = 0; i < nkeys; ++i)
  {
    perm |= dpf::gen(aeskey, point, dpfkey[i]) << i;
  }
  return perm;
}

template <uint8_t nservers>
inline __m256i partial_reduce(const __m256i & x)
{
  static const __m256i shuf0 = _mm256_set_epi8(
    (0xf0/nservers)*nservers, (0xe0/nservers)*nservers,
    (0xd0/nservers)*nservers, (0xc0/nservers)*nservers,
    (0xb0/nservers)*nservers, (0xa0/nservers)*nservers,
    (0x90/nservers)*nservers, (0x80/nservers)*nservers,
    (0x70/nservers)*nservers, (0x60/nservers)*nservers,
    (0x50/nservers)*nservers, (0x40/nservers)*nservers,
    (0x30/nservers)*nservers, (0x20/nservers)*nservers,
    (0x10/nservers)*nservers, (0x00/nservers)*nservers,
    (0xf0/nservers)*nservers, (0xe0/nservers)*nservers,
    (0xd0/nservers)*nservers, (0xc0/nservers)*nservers,
    (0xb0/nservers)*nservers, (0xa0/nservers)*nservers,
    (0x90/nservers)*nservers, (0x80/nservers)*nservers,
    (0x70/nservers)*nservers, (0x60/nservers)*nservers,
    (0x50/nservers)*nservers, (0x40/nservers)*nservers,
    (0x30/nservers)*nservers, (0x20/nservers)*nservers,
    (0x10/nservers)*nservers, (0x00/nservers)*nservers
  );

  return _mm256_sub_epi8(x, _mm256_shuffle_epi8(shuf0, _mm256_and_si256(_mm256_srli_epi16(x, 4), lo4)));
}

template <uint8_t nservers>
inline __m256i final_reduce(const __m256i & x)
{
  __m256i tmp;
  static constexpr uint8_t lg = static_cast<uint8_t>(std::ceil(std::log2(nservers)+1));
  if constexpr(lg > 5)
  {
    tmp = x;
    uint8_t * y = reinterpret_cast<uint8_t *>(&tmp);
    for (size_t i = 0; i < 32; ++i) y[i] %= nservers;
    return tmp;
  }
  else if constexpr(lg == 5)
  {
    static const __m256i y = _mm256_set1_epi8(static_cast<uint8_t>((0x10 / nservers) * nservers));
    __m256i tmp = _mm256_sub_epi8(x, _mm256_and_si256(_mm256_cmpgt_epi8(x, lo4), y));
  }
  else
  {
    tmp = x;
  }

  static const __m256i shuf = _mm256_set_epi8(
    (0x0f/nservers)*nservers, (0x0e/nservers)*nservers,
    (0x0d/nservers)*nservers, (0x0c/nservers)*nservers,
    (0x0b/nservers)*nservers, (0x0a/nservers)*nservers,
    (0x09/nservers)*nservers, (0x08/nservers)*nservers,
    (0x07/nservers)*nservers, (0x06/nservers)*nservers,
    (0x05/nservers)*nservers, (0x04/nservers)*nservers,
    (0x03/nservers)*nservers, (0x02/nservers)*nservers,
    (0x01/nservers)*nservers, (0x00/nservers)*nservers,
    (0x0f/nservers)*nservers, (0x0e/nservers)*nservers,
    (0x0d/nservers)*nservers, (0x0c/nservers)*nservers,
    (0x0b/nservers)*nservers, (0x0a/nservers)*nservers,
    (0x09/nservers)*nservers, (0x08/nservers)*nservers,
    (0x07/nservers)*nservers, (0x06/nservers)*nservers,
    (0x05/nservers)*nservers, (0x04/nservers)*nservers,
    (0x03/nservers)*nservers, (0x02/nservers)*nservers,
    (0x01/nservers)*nservers, (0x00/nservers)*nservers
  );
  
  return _mm256_sub_epi8(tmp, _mm256_shuffle_epi8(shuf, _mm256_and_si256(tmp, lo4)));
}

template <size_t nmasks, size_t nservers = static_cast<size_t>(std::exp2(nmasks)), bool do_reduce = static_cast<bool>(nservers < (std::exp2(nmasks)))>
inline static void splice(const __m128i mask[nmasks], __m1024i out)
{
  const static __m256i shuffle = _mm256_setr_epi64x(0x0000000000000000,
    0x0101010101010101, 0x0202020202020202, 0x0303030303030303);
  const static __m256i bit_mask = _mm256_set1_epi64x(0x7fbfdfeff7fbfdfe);
  const static __m256i ff = _mm256_set1_epi8(0xff);
  const __m256i shift[] = { _mm256_set1_epi8(1), _mm256_set1_epi8(2),
    _mm256_set1_epi8(4), _mm256_set1_epi8(8), _mm256_set1_epi8(16),
    _mm256_set1_epi8(32), _mm256_set1_epi8(64), _mm256_set1_epi8(128) };

  constexpr size_t lgbit = 1ULL << static_cast<size_t>(8 - std::ceil(std::log2(nservers)));

  for (size_t i = 0; i < nmasks; ++i)
  {
    const __m256i rem = (do_reduce) ? _mm256_set1_epi8(static_cast<uint8_t>((1ULL << i) % nservers)) : shift[i];
    const uint32_t * maski = (uint32_t *)&mask[i];
    for (int j = 0; j < 4; j++)
    {
      out[j] = _mm256_add_epi8(out[j], // accumulate into result
        _mm256_and_si256(rem,
          _mm256_cmpeq_epi8( // 0xff if byte is 0xff; else 0x00
            _mm256_or_si256( // set bits not possibly set by shuffle
              _mm256_shuffle_epi8( // shuffle 32 bits into 32 bytes
                _mm256_set1_epi32(maski[j]),
              shuffle),
            bit_mask),
          ff)
        )
      );
      //printf("before do_reduce!\n");
      if constexpr(do_reduce)
      {
        // printf("I've called!\n");
        if (i & lgbit) out[j] = partial_reduce<nservers>(out[j]);
      }
    }
  }
  if constexpr(do_reduce)
  {
    for (uint8_t j = 0; j < 4; ++j)
    {
      if constexpr(nmasks % lgbit) out[j] = partial_reduce<nservers>(out[j]);
      out[j] = final_reduce<nservers>(out[j]);
    }
  }
}

template <size_t nkeys, size_t nservers = static_cast<size_t>(std::exp2(nkeys)), size_t nitems>
inline void evalfull(AES_KEY & aeskey, dpf_key<nitems> * dpfkey,
  __m128i ** s, uint8_t ** t, uint8_t * output)
{
  constexpr size_t depth = dpf_key<nitems>::depth;
  constexpr size_t output_length = dpf_key<nitems>::output_length;

  __m128i child[2];
  uint8_t ts[2];
  for (size_t l = 0; l < nkeys; ++l)
  {
    int curlayer = depth % 2;

    __m128i * s_[2] = { s[l], s[l] + output_length/2 };
    uint8_t * t_[2] = { t[l], t[l] + output_length/2 };

    s_[curlayer][0] = dpfkey[l].root;
    t_[curlayer][0] = _mm_getlsb_si128(dpfkey[l].root);

    for (size_t i = 0; i < depth; ++i)
    {
      curlayer = 1 - curlayer;
      const size_t itemnumber = std::max(output_length >> (depth-i), 1UL);
      for (size_t j = 0; j < itemnumber; ++j)
      {
        expand(aeskey, s_[1-curlayer][j], child, ts);
        s_[curlayer][2*j] = _mm_xorif_si128(child[L], dpfkey[l].cw[i], t_[1-curlayer][j]);
        t_[curlayer][2*j] = ts[L] ^ dpfkey[l].t[i][L] & t_[1-curlayer][j];
        if (2*j+1 < 2*itemnumber)
        {
          s_[curlayer][2*j+1] = _mm_xorif_si128(child[R], dpfkey[l].cw[i], t_[1-curlayer][j]);
          t_[curlayer][2*j+1] = ts[R] ^ dpfkey[l].t[i][R] & t_[1-curlayer][j];
        }
      }
    }
  }

  memset(output, 0, nitems * sizeof(uint8_t));
  __m128i tmp[nkeys];
  __m1024i * output1024 = reinterpret_cast<__m1024i *>(output);
  for (size_t j = 0; j < output_length; ++j)
  {
    for (size_t l = 0; l < nkeys; ++l)
    {
      tmp[l] = _mm_xorif_si128(s[l][j], dpfkey[l].final, t[l][j]);
    }
    splice<nkeys,nservers>(tmp, output1024[j]);
  }
}

template <size_t nkeys, size_t nservers = static_cast<uint8_t>(std::exp2(nkeys)), size_t nitems>
inline size_t eval(AES_KEY & aeskey, dpf_key<nitems> dpfkey[nkeys],
  const size_t input)
{
  size_t result = 0;
  for (size_t i = 0; i < nkeys; ++i)
  {
  	result += (static_cast<size_t>(dpf::eval(aeskey, dpfkey[i], input)) << i) % nservers;
  }
  return result % nservers;
}

}

#endif