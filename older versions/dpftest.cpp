#include <cstdio>
#include <cstring>

#include "dpf.h"

int main(int argc, char * argv[])
{
  using namespace dpf;

  AES_KEY aeskey;
  AES_set_encrypt_key(_mm_set_epi64x(597349, 121379), &aeskey);

  typedef bool leaf_type;
  const size_t nitems = (1ULL << 12);
  const size_t item = 121;

  dpf::dpf_key<nitems, leaf_type> dpfkey[2];

  __m128i * s, * s1;
  if (posix_memalign((void**)&s, sizeof(__m128i), dpf_key<nitems,leaf_type>::output_length * sizeof(__m128i)) != 0)
    printf("alloc failed\n");
  if (posix_memalign((void**)&s1, sizeof(__m128i), dpf_key<nitems,leaf_type>::output_length * sizeof(__m128i)) != 0)
    printf("alloc failed\n");
  
  uint8_t * t = (uint8_t*)malloc(dpf_key<nitems,leaf_type>::output_length);
  uint8_t * t1 = (uint8_t*)malloc(dpf_key<nitems,leaf_type>::output_length);
  dpf::gen(aeskey, item, dpfkey, (leaf_type)1);
  
  dpf::evalfull(aeskey, dpfkey[0], s, t);
  dpf::evalfull(aeskey, dpfkey[1], s1, t1);
  printf("output_length: %d\n", (int)dpf_key<nitems,leaf_type>::output_length);
  printf("size of: %d\n", (int)sizeof(long long int));
  for(int i = 0; i < dpf_key<nitems,leaf_type>::output_length; ++i){
  	printf("%llx\n",  s[i][0] ^ s1[i][0]);
  	printf("%llx\n",  s[i][1] ^ s1[i][1]);
  }
   
  free(s);
  free(t);

  free(s1);
  free(t1);
  return 0;
  //1,0,0,0,0,0,0,0,0,0,0,0,0,0
}