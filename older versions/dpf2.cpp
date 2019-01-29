#include <algorithm>
#include <tuple>
#include <cmath>
#include "dpf2.h"

int main(int argc, char ** argv)
{
  AES_KEY aeskey;
  AES_set_encrypt_key(_mm_set_epi64x(597349, 121379), &aeskey);

  const size_t nitems = 1ULL << 10;
  const size_t item = 0;
  dpf_key<nitems> dpfkey[2];

  __m128i * s;
  if (posix_memalign((void**)&s, sizeof(__m128i), dpf_key<nitems>::length * sizeof(__m128i)) != 0)
  	printf("alloc failed\n");
  
  uint8_t * t = (uint8_t *)malloc(nitems * sizeof(uint8_t));

  gen(aeskey, item, dpfkey);
  
  evalfull(aeskey, dpfkey[0], s, t);
  for (int i = 0; i < dpf_key<nitems>::length; ++i) printf("%llx %llx ", s[i][0], s[i][1]); printf("\n");

  evalfull(aeskey, dpfkey[1], s, t);
  for (int i = 0; i < dpf_key<nitems>::length; ++i) printf("%llx %llx ", s[i][0], s[i][1]); printf("\n");
  
  free(s);
  free(t);

  return 0;
}