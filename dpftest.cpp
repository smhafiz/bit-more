#include <cstdio>
#include <cstring>

#include "dpf.h"

int main(int argc, char * argv[])
{
  using namespace dpf;

  AES_KEY aeskey;
  AES_set_encrypt_key(_mm_set_epi64x(597349, 121379), &aeskey);

  typedef bool leaf_type;
  const size_t nitems = (1ULL << 34);
  const size_t item = 100;

  dpf::dpf_key<nitems, leaf_type> dpfkey[2];

  __m128i * s;
  if (posix_memalign((void**)&s, sizeof(__m128i), dpf_key<nitems,leaf_type>::output_length * sizeof(__m128i)) != 0)
    printf("alloc failed\n");
  uint8_t * t = (uint8_t*)malloc(dpf_key<nitems,leaf_type>::output_length);

  dpf::gen(aeskey, item, dpfkey, (leaf_type)1);
  
  dpf::evalfull(aeskey, dpfkey[0], s, t);
  dpf::evalfull(aeskey, dpfkey[1], s, t);

  return 0;
}