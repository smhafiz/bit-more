#ifndef DPF_H__
#define DPF_H__

#include <cstdint>
#include <cmath>
#include <x86intrin.h>
#include <bsd/stdlib.h>
#include <cstring>

#include "aes.h"
#include "prg.h"

#define L 0
#define R 1
#define LG_128 7

static const __m128i lsb_mask[2] = {
  _mm_setzero_si128(),                                     // 0b00...0000
  _mm_set_epi64x(0,1),                                     // 0b00...0001
};
static const __m128i lsb_mask_inv = _mm_set_epi64x(-1,-2); // 0b11...1110
static const __m128i if_mask[2] = {
  _mm_setzero_si128(),                                     // 0b00...0000
  _mm_set1_epi8(-1)                                        // 0b11...1111
};

inline uint8_t _mm_getlsb_si128(const __m128i & x)
{
  __m128i vcmp = _mm_xor_si128(_mm_and_si128(x, lsb_mask[1]), lsb_mask[1]);
  return static_cast<uint8_t>(_mm_testz_si128(vcmp, vcmp));
}
inline __m128i _mm_clearlsb_si128(const __m128i & x)
{
  return _mm_and_si128(x, lsb_mask_inv);
}
inline __m128i _mm_setlsb_si128(const __m128i & x, bool b = true)
{
  return _mm_or_si128(_mm_clearlsb_si128(x), lsb_mask[b ? 1 : 0]);
}
inline __m128i _mm_xorif_si128(const __m128i & x, const __m128i & y, bool b)
{
  return _mm_xor_si128(x, _mm_and_si128(y, if_mask[b ? 1 : 0]));
}

struct cwbits { private: uint8_t t[2]; public: uint8_t & operator[](bool b) { return t[b ? 1 : 0]; } };
template<size_t nitems> struct dpf_key
{
  static constexpr size_t depth = static_cast<size_t>(std::ceil(std::log2(nitems / 128.0)));
  static constexpr size_t length = static_cast<size_t>(std::ceil(nitems / 128.0));

  __m128i root;
  __m128i cw[depth+1];
  cwbits t[depth];
};

inline void expand(AES_KEY & aeskey, const __m128i & seed, __m128i s[2],
  uint8_t t[2])
{
  const __m128i seedL = _mm_clearlsb_si128(seed);
  const __m128i seedR = _mm_setlsb_si128(seed);

  s[L] = seedL;
  s[R] = seedR;

  AES_ecb_encrypt_blks(s, 2, &aeskey);

  s[L] = _mm_xor_si128(s[L], seedL);
  t[L] = _mm_getlsb_si128(s[L]);

  s[R] = _mm_xor_si128(s[R], seedR);
  t[R] = _mm_getlsb_si128(s[R]);
}

inline void nextlayer(AES_KEY & aeskey, __m128i s[2], uint8_t t[2], uint8_t bit,
  __m128i s0[2], __m128i s1[2], uint8_t t0[2], uint8_t t1[2],
  __m128i & cw, cwbits & cwt, bool clearlsb = true)
{
//  __m128i s0[2], s1[2];
//  uint8_t t0[2], t1[2];

  expand(aeskey, s[0], s0, t0);
  expand(aeskey, s[1], s1, t1);

  const uint8_t keep = (bit == 0) ? L : R, lose = 1 - keep;
  cw = _mm_xor_si128(s0[lose], s1[lose]);
  if (clearlsb) cw = _mm_clearlsb_si128(cw);

  cwt[L] = t0[L] ^ t1[L] ^ bit ^ 1;
  cwt[R] = t0[R] ^ t1[R] ^ bit;

  s[L] = _mm_xorif_si128(s0[keep], cw, t[L]);
  t[L] = t0[keep] ^ (t[L] & cwt[keep]);
  s[R] = _mm_xorif_si128(s1[keep], cw, t[R]);
  t[R] = t1[keep] ^ (t[R] & cwt[keep]);
}

template <size_t nitems>
void gen(AES_KEY & aeskey, size_t point, dpf_key<nitems> dpfkey[2])
{
  __m128i s[2], s0[2], s1[2];
  uint8_t t[2], t0[2], t1[2];
  arc4random_buf(s, 2 * sizeof(__m128i));
  t[0] = _mm_getlsb_si128(s[0]);
  dpfkey[0].root = s[0];
  t[1] = !t[0];
  dpfkey[1].root = _mm_setlsb_si128(s[1], t[1]);

  for (size_t i = 0; i < dpf_key<nitems>::depth - 1; ++i)
  {
    const uint8_t bit = (point >> (dpf_key<nitems>::depth - i - 1)) & 1U;
    nextlayer(aeskey, s, t, bit, s0, s1, t0, t1, dpfkey[0].cw[i], dpfkey[0].t[i]);
  }
  const uint8_t bit = (point >> LG_128) & 1U;
  const uint8_t keep = (bit == 0) ? L : R;
  const uint8_t lose = 1 - keep;
  const size_t lo = (point & 0x40) ? 0ULL : 1ULL;
  
  nextlayer(aeskey, s, t, bit, s0, s1, t0, t1, dpfkey[0].cw[dpf_key<nitems>::depth-1],
    dpfkey[0].t[dpf_key<nitems>::depth-1], false);

  __m128i finalmask = _mm_slli_epi64(_mm_set_epi64x(1-lo, lo), point % 0x40);
  dpfkey[0].cw[dpf_key<nitems>::depth] = _mm_xor_si128(finalmask, _mm_xor_si128(
    _mm_xorif_si128(s0[keep], dpfkey[0].cw[dpf_key<nitems>::depth-1], t[L]),
    _mm_xorif_si128(s1[keep], dpfkey[0].cw[dpf_key<nitems>::depth-1], t[R]))
  );

  memcpy(&dpfkey[1].cw, &dpfkey[0].cw, sizeof(dpf_key<nitems>::cw));
  memcpy(&dpfkey[1].t, &dpfkey[0].t, sizeof(dpf_key<nitems>::t));
}

template <size_t nitems>
inline void dolayer(AES_KEY & aeskey,__m128i * s, uint8_t * t,
  const size_t stepsize, dpf_key<nitems> dpfkey, size_t i,
  __m128i child[2], uint8_t ts[2])
{
  for (size_t j = 0; j < dpf_key<nitems>::length; j += 2*stepsize)
  {
    expand(aeskey, s[j], child, ts);
    if (j+stepsize < nitems)
    {
      s[j+stepsize] = _mm_xorif_si128(child[R], dpfkey.cw[i], t[j]);
      t[j+stepsize] = ts[R] ^ dpfkey.t[i][R] & t[j];
    }
    s[j] = _mm_xorif_si128(child[L], dpfkey.cw[i], t[j]);
    t[j] = ts[L] ^ dpfkey.t[i][L] & t[j];
  }
}

template <size_t nitems>
inline void evalfull(AES_KEY & aeskey, dpf_key<nitems> & dpfkey, __m128i * s,
  uint8_t * t)
{
  s[0] = dpfkey.root;
  t[0] = _mm_getlsb_si128(dpfkey.root);

  __m128i child[2];
  uint8_t ts[2];
  size_t stepsize = 1ULL << (dpf_key<nitems>::depth - 1);
  for (size_t i = 0; i < dpf_key<nitems>::depth; ++i, stepsize /= 2)
  {
    dolayer(aeskey, s, t, stepsize, dpfkey, i, child, ts);
  }

  for (size_t j = 0; j < dpf_key<nitems>::length; ++j)
  {
    s[j] = _mm_xorif_si128(s[j], dpfkey.cw[dpf_key<nitems>::depth], t[j]);
  }
}

#endif