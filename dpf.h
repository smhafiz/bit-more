// evalrange with in place and alternating eval
//
// mux (pir) and demux; determine counts efficiently
//
// proof of consistency via batch testing, power trick

#ifndef DPF_H__
#define DPF_H__

#include <bsd/stdlib.h> // for arc4random_buf
#include <stdexcept>    // for std::runtime_error
#include <cmath>        // for std::log2 and std::ceil
#include <climits>      // for CHAR_BIT
#include <cstring>      // for std::memcpy

#include <x86intrin.h>

#include "prg.h"
#include "aes.h"

#define L 0
#define R 1

namespace dpf
{

static const __m128i lsb_mask[2] = {
  _mm_setzero_si128(),                                     // 0b00...0000
  _mm_set_epi64x(0,1),                                     // 0b00...0001
};
static const __m128i lsb_mask_inv = _mm_set_epi64x(-1,-2); // 0b11...1100
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

template<size_t nitems, typename leaf_type = bool, bool B = (sizeof(leaf_type) < 16)>
struct dpf_key;

struct cwbits
{
private:
  uint8_t t[2];
public:
  inline uint8_t & operator[](bool b) { return t[b ? 1U : 0U]; }
};

template<size_t nitems, typename leaf_type>
struct dpf_key<nitems, leaf_type, true>
{
  static constexpr size_t leaf_bits = std::is_same<leaf_type, bool>::value ? 1 : sizeof(leaf_type) * 8;
  static constexpr size_t input_bits = static_cast<size_t>(std::ceil(std::log2(nitems)));
  static constexpr size_t outs_per_m128 = static_cast<size_t>(128.0 / leaf_bits);
  static constexpr size_t depth = static_cast<size_t>(std::ceil(std::log2(nitems / outs_per_m128)));
  static constexpr size_t output_length = static_cast<size_t>(std::ceil(nitems * leaf_bits / 128.0));
  static constexpr size_t lg_outs_per_m128 = static_cast<size_t>(std::log2(outs_per_m128));

  static_assert(nitems > outs_per_m128);
  static_assert(128 % leaf_bits == 0);

  __m128i root;
  __m128i cw[depth];
  __m128i final;
  cwbits t[depth];  
};

template<size_t nitems, typename leaf_type>
struct dpf_key<nitems, leaf_type, false>
{
  static constexpr size_t leaf_bits = std::is_same<leaf_type, bool>::value ? 1 : CHAR_BIT * sizeof(leaf_type);
  static_assert(leaf_bits % 128 == 0);
  static constexpr size_t input_bits = static_cast<size_t>(std::ceil(std::log2(nitems)));
  static constexpr size_t depth = static_cast<size_t>(std::ceil(std::log2(nitems / leaf_bits)));
  static constexpr size_t output_bytes = static_cast<size_t>(std::ceil(nitems * leaf_bits / 8.0));
  static constexpr size_t m128s_per_out = static_cast<size_t>(leaf_bits / 128);

  __m128i root;
  __m128i cw[depth];
  __m128i final[m128s_per_out];
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
  expand(aeskey, s[0], s0, t0);
  expand(aeskey, s[1], s1, t1);

  const uint8_t keep = (bit == 0) ? L : R, lose = 1 - keep;
  cw = _mm_xor_si128(s0[lose], s1[lose]);
  if (clearlsb) cw = _mm_clearlsb_si128(cw);

  cwt[L] = t0[L] ^ t1[L] ^ !bit;
  cwt[R] = t0[R] ^ t1[R] ^ bit;

  s[L] = _mm_xorif_si128(s0[keep], cw, t[L]);
  t[L] = t0[keep] ^ (t[L] & cwt[keep]);
  s[R] = _mm_xorif_si128(s1[keep], cw, t[R]);
  t[R] = t1[keep] ^ (t[R] & cwt[keep]);
}

template <size_t nitems, typename leaf_type>
void doleaf(size_t point, leaf_type output, dpf_key<nitems,leaf_type> & dpfkey,
  __m128i s0[2], __m128i s1[2], uint8_t t[2])
{
  constexpr size_t depth = dpf_key<nitems,leaf_type>::depth;
  constexpr size_t lg_outs_per_m128 = dpf_key<nitems,leaf_type>::lg_outs_per_m128;
  constexpr size_t outs_per_u64 = dpf_key<nitems,leaf_type>::outs_per_m128 / 2;
  constexpr size_t leaf_bits = dpf_key<nitems,leaf_type>::leaf_bits;

  const uint8_t keep = (((point >> lg_outs_per_m128) & 1U) == 0) ? L : R;
  const size_t lo = (point & outs_per_u64) ? 0ULL : 1ULL;
  dpfkey.final = _mm_set_epi64x(output*(1-lo), output*lo);
  dpfkey.final = _mm_slli_epi64(dpfkey.final, leaf_bits * (point % outs_per_u64));
  dpfkey.final = _mm_xor_si128(dpfkey.final, _mm_xor_si128(
    _mm_xorif_si128(s0[keep], dpfkey.cw[depth-1], t[L]),
    _mm_xorif_si128(s1[keep], dpfkey.cw[depth-1], t[R]))
  );
}

template <size_t nitems, typename leaf_type>
void gen(AES_KEY & aeskey, size_t point, dpf_key<nitems,leaf_type> dpfkey[2],
  leaf_type output)
{
  if (point >= nitems)
  {
    throw std::runtime_error("point is out of range");
  }
  constexpr size_t depth = dpf_key<nitems,leaf_type>::depth;
  constexpr size_t input_bits = dpf_key<nitems,leaf_type>::input_bits;

  __m128i s[2], s0[2], s1[2];
  uint8_t t[2], t0[2], t1[2];
  arc4random_buf(s, 2 * sizeof(__m128i));
  t[0] = _mm_getlsb_si128(s[0]);
  dpfkey[0].root = s[0];
  t[1] = !t[0];
  dpfkey[1].root = _mm_setlsb_si128(s[1], t[1]);

  for (size_t i = 0; i < depth; ++i)
  {
    const uint8_t bit = (point >> (input_bits-i-1)) & 1U;
    nextlayer(aeskey, s, t, bit, s0, s1, t0, t1, dpfkey[0].cw[i], dpfkey[0].t[i],
      (i < depth - 1));
  }

  doleaf(point, output, dpfkey[0], s0, s1, t);

  memcpy(&dpfkey[1].cw, &dpfkey[0].cw, sizeof(dpf_key<nitems,leaf_type>::cw));
  memcpy(&dpfkey[1].final, &dpfkey[0].final, sizeof(dpf_key<nitems,leaf_type>::final));
  memcpy(&dpfkey[1].t, &dpfkey[0].t, sizeof(dpf_key<nitems,leaf_type>::t));
}

/*template <size_t nitems, typename leaf_type>
inline void evalfull2(AES_KEY & aeskey, dpf_key<nitems,leaf_type> & dpfkey,
  __m128i * s, uint8_t * t)
{
  constexpr size_t depth = dpf_key<nitems,leaf_type>::depth;
  constexpr size_t output_length = dpf_key<nitems,leaf_type>::output_length;

  s[0] = dpfkey.root;
  t[0] = _mm_getlsb_si128(dpfkey.root);

  size_t stepsize = 1ULL << (depth - 1);
  for (size_t i = 0; i < depth; ++i, stepsize /= 2)
  {
    uint8_t ts[2];
    for (size_t j = 0; j < output_length; j += 2*stepsize)
    {
      expand(aeskey, s[j], &s[j], ts);
      if (j+stepsize < output_length)
      {
        s[j+stepsize] = _mm_xorif_si128(s[j+1], dpfkey.cw[i], t[j]);
        t[j+stepsize] = ts[R] ^ (dpfkey.t[i][R] & t[j]);
      }
      s[j] = _mm_xorif_si128(s[j], dpfkey.cw[i], t[j]);
      t[j] = ts[L] ^ (dpfkey.t[i][L] & t[j]);
    }
  }

  for (size_t j = 0; j < output_length; ++j)
  {
    s[j] = _mm_xorif_si128(s[j], dpfkey.final, t[j]);
  }
}*/

template <size_t nitems, typename leaf_type>
inline void evalfull(AES_KEY & aeskey, dpf_key<nitems,leaf_type> & dpfkey,
  __m128i * s, uint8_t * t)
{
  constexpr size_t depth = dpf_key<nitems,leaf_type>::depth;
  constexpr size_t output_length = dpf_key<nitems,leaf_type>::output_length;

  __m128i * s_[2] = { s, s + output_length/2 };
  uint8_t * t_[2] = { t, t + output_length/2 };

  int curlayer = depth % 2;

  s_[curlayer][0] = dpfkey.root;
  t_[curlayer][0] = _mm_getlsb_si128(dpfkey.root);

  __m128i child[2];
  uint8_t ts[2];
  for (size_t i = 0; i < depth; ++i)
  {
    curlayer = 1 - curlayer;
    const size_t itemnumber = std::max((output_length / (1ULL << (depth - i))), 1ULL);
    for (size_t j = 0; j < itemnumber; ++j)
    {
      expand(aeskey, s_[1-curlayer][j], child, ts);
      s_[curlayer][2*j] = _mm_xorif_si128(child[L], dpfkey.cw[i], t_[1-curlayer][j]);
      t_[curlayer][2*j] = ts[L] ^ dpfkey.t[i][L] & t_[1-curlayer][j];
      if (2*j+1 < 2*itemnumber)
      {
        s_[curlayer][2*j+1] = _mm_xorif_si128(child[R], dpfkey.cw[i], t_[1-curlayer][j]);
        t_[curlayer][2*j+1] = ts[R] ^ dpfkey.t[i][R] & t_[1-curlayer][j];
      }
    }
  }

  for (size_t j = 0; j < output_length; ++j)
  {
    s[j] = _mm_xorif_si128(s[j], dpfkey.final, t[j]);
  }
}

template<size_t nitems>
inline bool getword(dpf_key<nitems, bool> dpfkey, const size_t input,
  __m128i & S, uint8_t T)
{
  const size_t lo = (input & 64) ? 0ULL : 1ULL;
  const __m128i mask = _mm_slli_epi64(_mm_set_epi64x(1-lo, lo), input % 64);
  S = _mm_and_si128(_mm_xorif_si128(S, dpfkey.final, T), mask);

  return static_cast<bool>(_mm_testz_si128(S, mask));
}

template<size_t nitems, typename leaf_type>
inline leaf_type getword(dpf_key<nitems,leaf_type,true> dpfkey, const size_t input,
  __m128i & S, uint8_t T)
{
  constexpr size_t outs_per_u64 = dpf_key<nitems,leaf_type>::outs_per_m128 / 2;
  constexpr size_t leaf_bits = dpf_key<nitems,leaf_type>::leaf_bits;

  S = _mm_srli_epi64(_mm_xorif_si128(S, dpfkey.final, T), leaf_bits*(input % outs_per_u64));
  if constexpr(sizeof(leaf_type) == 1)
  {
    return static_cast<leaf_type>(input & outs_per_u64 ? _mm_extract_epi8(S, outs_per_u64) : _mm_extract_epi8(S, 0));
  }
  else if constexpr(sizeof(leaf_type) == 2)
  {
    return static_cast<leaf_type>(input & outs_per_u64 ? _mm_extract_epi16(S, outs_per_u64) : _mm_extract_epi16(S, 0));
  }
  else if constexpr(sizeof(leaf_type) == 4)
  {
    return static_cast<leaf_type>(input & outs_per_u64 ? _mm_extract_epi32(S, outs_per_u64) : _mm_extract_epi32(S, 0));
  }
  else if constexpr(sizeof(leaf_type) == 8)
  {
    return static_cast<leaf_type>(input & outs_per_u64 ? _mm_extract_epi64(S, outs_per_u64) : _mm_extract_epi64(S, 0));
  }
  else if constexpr(sizeof(leaf_type) == 16)
  {
    return static_cast<leaf_type>(S);
  }
}

template <size_t nitems, typename leaf_type>
inline leaf_type eval(AES_KEY & aeskey, dpf_key<nitems, leaf_type> & dpfkey,
  const size_t input)
{
  constexpr size_t depth = dpf_key<nitems,leaf_type>::depth;
  constexpr size_t input_bits = dpf_key<nitems,leaf_type>::input_bits;

  __m128i S = dpfkey.root;
  uint8_t T = _mm_getlsb_si128(S);
  __m128i child[2];
  uint8_t ts[2];
  for (size_t i = 0; i < depth; ++i)
  {
    const uint8_t bit = (input >> (input_bits-i-1)) & 1U;
    expand(aeskey, S, child, ts); 

    S = _mm_xorif_si128(child[bit], dpfkey.cw[i], T);
    T = ts[bit] ^ (dpfkey.t[i][bit] & T);
  }

  return getword(dpfkey, input, S, T);
}

};

#endif