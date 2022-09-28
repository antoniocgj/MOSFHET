#include "mosfhet.h"

// key expansion, adapted from: https://www.intel.com/content/dam/doc/white-paper/advanced-encryption-standard-new-instructions-set-paper.pdf
static inline __m128i AES_128_ASSIST (__m128i temp1, __m128i temp2){
  __m128i temp3;
  temp2 = _mm_shuffle_epi32 (temp2 ,0xff);
  temp3 = _mm_slli_si128 (temp1, 0x4);
  temp1 = _mm_xor_si128 (temp1, temp3);
  temp3 = _mm_slli_si128 (temp3, 0x4);
  temp1 = _mm_xor_si128 (temp1, temp3);
  temp3 = _mm_slli_si128 (temp3, 0x4);
  temp1 = _mm_xor_si128 (temp1, temp3);
  temp1 = _mm_xor_si128 (temp1, temp2);
  return temp1;
}

#ifdef USE_VAES
static __m512i __global_prng_aes_key[11];
#else
static __m128i __global_prng_aes_key[11];
#define _mm512_broadcast_i64x2(X) X 
#endif


void setup_aes_prgn_key(__m128i * seed){
  const uint8_t rnd_c[11] = {0, 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1B, 0x36};
  __m128i tmp1, tmp2;
  tmp1 = _mm_loadu_si128((__m128i*)seed);
  __global_prng_aes_key[0] = _mm512_broadcast_i64x2(tmp1);
  for (size_t i = 1; i < 11; i++){
    tmp2 = _mm_aeskeygenassist_si128(tmp1, rnd_c[i]);
    tmp1 = AES_128_ASSIST(tmp1, tmp2);
    __global_prng_aes_key[i] = _mm512_broadcast_i64x2(tmp1);
  }
}

void aes_prgn_setup_rnd_seed(){
  __m128i s[2];
  generate_rnd_seed((uint64_t*)s);
  setup_aes_prgn_key(s);
}

#ifdef USE_VAES
void aes_prgn_next_16(__m512i * out, __m512i * cnt){
  const __m512i incv = {0, 4, 0, 4, 0, 4, 0, 4};
  for (size_t i = 0; i < 4; i++){
    out[i] = _mm512_xor_si512(*cnt, __global_prng_aes_key[0]);
    *cnt = _mm512_add_epi64 (*cnt, incv);
  }

  for (size_t i = 1; i < 10; i++){
    for (size_t j = 0; j < 4; j++){
      out[j] = _mm512_aesenc_epi128(out[j], __global_prng_aes_key[i]);
    }
  }

  for (size_t i = 0; i < 4; i++){
    out[i] = _mm512_aesenclast_epi128(out[i], __global_prng_aes_key[10]);
  }
}

// Generates outlen bytes from the first 8 bytes of *input
// Assumes: output is aligned, outlen >= 256 , inlen >= 16
void aes_prng(uint8_t *output, uint64_t outlen, const uint8_t *input,  uint64_t inlen){
  assert(outlen >= 256);
  assert(inlen >= 16);
  static bool first_exec = true;
  if(first_exec) {
    aes_prgn_setup_rnd_seed();
    first_exec = false;
  }
  __m128i cnt = _mm_loadu_si128((__m128i *)input);
  __m512i cntv = _mm512_broadcast_i64x2 (cnt);
  __m512i cntv2 = {0, 0, 0, 1, 0, 2, 0, 3};
  cntv = _mm512_add_epi64 (cntv, cntv2);
  size_t i;
  for (i = 0; i < outlen - 255; i+=256){
    aes_prgn_next_16((__m512i *) &(output[i]), &cntv);
  }
  if(outlen > i){
    __m512i tmp[4];
    aes_prgn_next_16(tmp, &cntv);
    memcpy(&output[i], tmp, outlen - i);
  }
}

#else
void aes_prgn_next_16(__m128i * out, __m128i * cnt){
  const __m128i incv = {0, 4};
  for (size_t i = 0; i < 16; i++){
    out[i] = _mm_xor_si128(*cnt, __global_prng_aes_key[0]);
    *cnt = _mm_add_epi64 (*cnt, incv);
  }

  for (size_t i = 1; i < 10; i++){
    for (size_t j = 0; j < 16; j++){
      out[j] = _mm_aesenc_si128(out[j], __global_prng_aes_key[i]);
    }
  }

  for (size_t i = 0; i < 16; i++){
    out[i] = _mm_aesenclast_si128 (out[i], __global_prng_aes_key[10]);
  }
}

// Generates outlen bytes from the first 8 bytes of *input
// Assumes: output is aligned, outlen >= 256 , inlen >= 16
void aes_prng(uint8_t *output, uint64_t outlen, const uint8_t *input,  uint64_t inlen){
  assert(outlen >= 256);
  assert(inlen >= 16);
  static bool first_exec = true;
  if(first_exec) {
    aes_prgn_setup_rnd_seed();
    first_exec = false;
  }
  __m128i cnt = _mm_loadu_si128((__m128i *)input);
  size_t i;
  for (i = 0; i < outlen - 255; i+=256){
    aes_prgn_next_16((__m128i *) &(output[i]), &cnt);
  }
  if(outlen > i){
    __m128i tmp[16];
    aes_prgn_next_16(tmp, &cnt);
    memcpy(&output[i], tmp, outlen - i);
  }
}
#endif

