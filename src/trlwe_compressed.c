#include "mosfhet.h"
// #define USE_SHAKE
// Functions from xoroshiro128++ 1.0, by David Blackman and Sebastiano Vigna (vigna@acm.org)
// Xoroshiro License: See <http://creativecommons.org/publicdomain/zero/1.0/>. 
inline uint64_t rotl(const uint64_t x, int k) {
	return (x << k) | (x >> (64 - k));
}

uint64_t xoroshiro128pp_next(uint64_t s[2]) {
	const uint64_t s0 = s[0];
	uint64_t s1 = s[1];
	const uint64_t result = rotl(s0 + s1, 17) + s0;

	s1 ^= s0;
	s[0] = rotl(s0, 49) ^ s1 ^ (s1 << 21); // a, b
	s[1] = rotl(s1, 28); // c

	return result;
}

void xoroshiro128pp_vnext(uint64_t * const restrict array, uint64_t sv[2][4]) {
 	uint64_t s1[4];
  for(int i = 0; i < 4; i++) array[i] = rotl(sv[0][i] + sv[1][i], 17) + sv[0][i];
  for(int i = 0; i < 4; i++) s1[i] = sv[0][i] ^ sv[1][i];
  for(int i = 0; i < 4; i++) sv[0][i] = rotl(sv[0][i], 49) ^ s1[i] ^ (s1[i] << 21);
  for(int i = 0; i < 4; i++) sv[1][i] = rotl(s1[i], 28);
}

inline void xoroshiro128pp_vnext_sub(uint64_t * const restrict array, uint64_t sv[2][4]) {
 	uint64_t s1[4];
  for(int i = 0; i < 4; i++) array[i] -= rotl(sv[0][i] + sv[1][i], 17) + sv[0][i];
  for(int i = 0; i < 4; i++) s1[i] = sv[0][i] ^ sv[1][i];
  for(int i = 0; i < 4; i++) sv[0][i] = rotl(sv[0][i], 49) ^ s1[i] ^ (s1[i] << 21);
  for(int i = 0; i < 4; i++) sv[1][i] = rotl(s1[i], 28);
}

// Shake128
void shake128_absorb(uint64_t *s, const unsigned char *input, unsigned int inputByteLen);
void shake128(unsigned char *output, unsigned long long outlen, const unsigned char *input,  unsigned long long inlen);
void shake128_squeeze(unsigned char *output, unsigned long long outlen, uint64_t * s);

// TRLWE Compressed Sample Functions
// 16-byte random seed;

TRLWE trlwe_alloc_new_compressed_sample(int k, int N){
  TRLWE res;
  res = (TRLWE) safe_malloc(sizeof(*res));
  res->a = (TorusPolynomial *) safe_malloc(sizeof(TorusPolynomial));
  res->a[0] = polynomial_new_torus_polynomial(16/sizeof(Torus));
  res->b = polynomial_new_torus_polynomial(N);
  res->k = k;
  return res;
}

void trlwe_load_compressed_sample(FILE * fd, TRLWE c){
  const int N = c->b->N;
  fread(c->a[0]->coeffs, sizeof(Torus), 16, fd);
  fread(c->b->coeffs, sizeof(Torus), N, fd);
}

TRLWE trlwe_load_new_compressed_sample(FILE * fd, int k, int N){
  TRLWE res = trlwe_alloc_new_compressed_sample(k, N);
  trlwe_load_compressed_sample(fd, res);
  return res;
}

void trlwe_save_compressed_sample(FILE * fd, TRLWE c){
  fwrite(c->a[0]->coeffs, sizeof(Torus), 16, fd);
  fwrite(c->b->coeffs, sizeof(Torus), c->b->N, fd);
}

void trlwe_compressed_sample(TRLWE out, TorusPolynomial m, TRLWE_Key key){
  const int N = key->s[0]->N, byte_size = 16;
  TorusPolynomial p_tmp = polynomial_new_torus_polynomial(N);

  generate_random_bytes(byte_size, (uint8_t *) out->a[0]->coeffs);
  uint8_t seedi[byte_size];
  memcpy(seedi, out->a[0]->coeffs, byte_size);
#ifdef USE_SHAKE
  uint64_t s[25] = {0};
  shake128_absorb(s, seedi, 16);
#else
  uint64_t seed[2][4];
  for (size_t i = 0; i < 4; i++){
    seed[0][i] = xoroshiro128pp_next((uint64_t *) seedi);
    seed[1][i] = xoroshiro128pp_next((uint64_t *) seedi);
  }
#endif

  // add error
  generate_torus_normal_random_array(out->b->coeffs, key->sigma, N);

  // internal product
  for (size_t i = 0; i < key->k; i++){
  #ifdef USE_SHAKE
    shake128_squeeze((uint8_t *) p_tmp->coeffs, sizeof(Torus)*N, s);
  #else
    for (size_t j = 0; j < N; j+=4) xoroshiro128pp_vnext(&p_tmp->coeffs[j], seed);
  #endif
    polynomial_naive_mul_addto_torus(out->b, p_tmp, key->s[i]);
  }

  if(m != NULL){
    for (size_t i = 0; i < N; i++){
      out->b->coeffs[i] += m->coeffs[i];
    }
  }
  free_polynomial(p_tmp);
}

TRLWE trlwe_new_compressed_sample(TorusPolynomial m, TRLWE_Key key){
  const int N = key->s[0]->N;
  TRLWE res = trlwe_alloc_new_compressed_sample(key->k, N);
  trlwe_compressed_sample(res, m, key);
  return res;
}

#if defined(AVX2_OPT) && !defined(USE_SHAKE)
static inline __m256i rotl_v(const __m256i x, const int k) {
	return _mm256_or_si256(_mm256_slli_epi64(x, k), _mm256_srli_epi64(x, 64 - k));
}

void trlwe_compressed_subto(TRLWE out, TRLWE in){
  const int N = out->b->N, k = out->k;
  uint8_t seedi[16];
  memcpy(seedi, in->a[0]->coeffs, 16);
  __m256i seed0, seed1;
  for (size_t i = 0; i < 4; i++){
    ((uint64_t * )&seed0)[i] = xoroshiro128pp_next((uint64_t *) seedi);
    ((uint64_t * )&seed1)[i] = xoroshiro128pp_next((uint64_t *) seedi);
  }
  __m256i s1;
  __m256i * outav = (__m256i *) out->a[0]->coeffs;
  __m256i * outbv = (__m256i *) out->b->coeffs;
  __m256i * inbv = (__m256i *) in->b->coeffs;

  for (size_t j = 0; j < N/4; j++){
    outav[j] -= rotl_v(seed0 + seed1, 17) + seed0;
    s1 = seed0 ^ seed1;
    seed0 = rotl_v(seed0, 49) ^ s1 ^ (s1 << 21);
    seed1 = rotl_v(s1, 28);
    outbv[j] -= inbv[j]; 
  }

  for (size_t i = 1; i < k; i++){ // Optimized for k = 1.
    outav = (__m256i *) out->a[k]->coeffs;
    for (size_t j = 0; j < N/4; j++){
      outav[j] -= rotl_v(seed0 + seed1, 17) + seed0;
      s1 = seed0 ^ seed1;
      seed0 = rotl_v(seed0, 49) ^ s1 ^ (s1 << 21);
      seed1 = rotl_v(s1, 28);
    }
  }
}

#else

void trlwe_compressed_subto(TRLWE out, TRLWE in){
  const int N = out->b->N, k = out->k;
  uint8_t seedi[16];
  memcpy(seedi, in->a[0]->coeffs, 16);

#ifndef USE_SHAKE
  uint64_t seed[2][4];
  for (size_t i = 0; i < 4; i++){
    seed[0][i] = xoroshiro128pp_next((uint64_t *) seedi);
    seed[1][i] = xoroshiro128pp_next((uint64_t *) seedi);
  }

  for (size_t j = 0; j < N; j+=4){
    xoroshiro128pp_vnext_sub(&out->a[0]->coeffs[j], seed);
    for (size_t i = 0; i < 4; i++) out->b->coeffs[j + i] -= in->b->coeffs[j + i]; 
  }
  for (size_t i = 1; i < k; i++){
    for (size_t j = 0; j < N; j+=4) xoroshiro128pp_vnext_sub(&out->a[i]->coeffs[j], seed);
  }
#else
  uint64_t s[25] = {0};
  shake128_absorb(s, seedi, 16);
  Torus tmp[N];
  for (size_t i = 0; i < k; i++){
    shake128_squeeze((uint8_t *) tmp, sizeof(Torus)*N, s);
    for (size_t j = 0; j < N; j++) out->a[i]->coeffs[j] -= tmp[j];
  }
  for (size_t j = 0; j < N; j++) out->b->coeffs[j] -= in->b->coeffs[j];
#endif

}

#endif



