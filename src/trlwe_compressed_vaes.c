#include "mosfhet.h"
// k = 1
// aes setup
#define ID_SIZE 8
// key expansion, adapted from: https://www.intel.com/content/dam/doc/white-paper/advanced-encryption-standard-new-instructions-set-paper.pdf
inline __m128i AES_128_ASSIST (__m128i temp1, __m128i temp2){
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

__m512i __global_aes_key[11];
void setup_aes_key(__m128i * seed){
  const uint8_t rnd_c[11] = {0, 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1B, 0x36};
  __m128i tmp1, tmp2;
  tmp1 = _mm_loadu_si128((__m128i*)seed);
  __global_aes_key[0] = _mm512_broadcast_i64x2(tmp1);
  for (size_t i = 1; i < 11; i++){
    tmp2 = _mm_aeskeygenassist_si128(tmp1, rnd_c[i]);
    tmp1 = AES_128_ASSIST(tmp1, tmp2);
    __global_aes_key[i] = _mm512_broadcast_i64x2(tmp1);
  }
}

void aes_setup_rnd_seed(){
  __m128i s;
  generate_random_bytes(sizeof(__m128i), (uint8_t*)&s);
  setup_aes_key(&s);
}

void aes_next_16(__m512i * out, __m512i * cnt){
  const __m512i incv = {0, 4, 0, 4, 0, 4, 0, 4};
  for (size_t i = 0; i < 4; i++){
    out[i] = _mm512_xor_si512(*cnt, __global_aes_key[0]);
    *cnt = _mm512_add_epi64 (*cnt, incv);
  }

  for (size_t i = 1; i < 10; i++){
    for (size_t j = 0; j < 4; j++){
      out[j] = _mm512_aesenc_epi128(out[j], __global_aes_key[i]);
    }
  }

  for (size_t i = 0; i < 4; i++){
    out[i] = _mm512_aesenclast_epi128(out[i], __global_aes_key[10]);
  }
}

void aes_next_4(__m512i * out, __m512i * cnt){
  const __m512i incv = {0, 4, 0, 4, 0, 4, 0, 4};
  *out = _mm512_xor_si512(*cnt, __global_aes_key[0]);
  *cnt = _mm512_add_epi64 (*cnt, incv);

  for (size_t i = 1; i < 10; i++){
    *out = _mm512_aesenc_epi128(*out, __global_aes_key[i]);
  }

  *out = _mm512_aesenclast_epi128(*out, __global_aes_key[10]);
}

// size is the number of 64-bit words to generate
void aes_full_expand(__m512i * out, uint64_t id, uint64_t size){
  size /= (512/64); // number of 512-bit words
  __m512i idv =  _mm512_maskz_set1_epi64 (0x55, id);
  __m512i cntv = {0, 0, 0, 1, 0, 2, 0, 3};
  cntv = _mm512_add_epi64 (cntv, idv);
  for (size_t i = 0; i < size; i+=4){
    aes_next_16(&out[i], &cntv);
  }
}


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
  fread(c->a[0]->coeffs, sizeof(Torus), ID_SIZE, fd);
  fread(c->b->coeffs, sizeof(Torus), N, fd);
}

TRLWE trlwe_load_new_compressed_sample(FILE * fd, int k, int N){
  TRLWE res = trlwe_alloc_new_compressed_sample(k, N);
  trlwe_load_compressed_sample(fd, res);
  return res;
}

void trlwe_save_compressed_sample(FILE * fd, TRLWE c){
  fwrite(c->a[0]->coeffs, sizeof(Torus), ID_SIZE, fd);
  fwrite(c->b->coeffs, sizeof(Torus), c->b->N, fd);
}

void trlwe_compressed_sample(TRLWE out, TorusPolynomial m, TRLWE_Key key){
  const int N = key->s[0]->N, byte_size = ID_SIZE;
  TorusPolynomial p_tmp = polynomial_new_torus_polynomial(N);

  generate_random_bytes(byte_size, (uint8_t *) out->a[0]->coeffs);
  uint64_t seedi = out->a[0]->coeffs[0];

  // add error
  generate_torus_normal_random_array(out->b->coeffs, key->sigma, N);

  // internal product
  for (size_t i = 0; i < key->k; i++){
    aes_full_expand((__m512i *)p_tmp->coeffs, seedi, N);
    polynomial_mul_addto_torus(out->b, p_tmp, key->s[i]);
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


void trlwe_compressed_subto(TRLWE out, TRLWE in){
  const int N = out->b->N;

  const uint64_t id = in->a[0]->coeffs[0];
  __m512i * outav = (__m512i *) out->a[0]->coeffs;
  __m512i * outbv = (__m512i *) out->b->coeffs;
  __m512i * inbv = (__m512i *) in->b->coeffs;
  __m512i idv =  _mm512_maskz_set1_epi64 (0x55, id);
  __m512i cntv = {0, 0, 0, 1, 0, 2, 0, 3};
  cntv = _mm512_add_epi64 (cntv, idv);
  __m512i tmp[4];

  for (size_t j = 0; j < N/8; j+=4){
    aes_next_16(tmp, &cntv);
    for (size_t i = 0; i < 4; i++){
      outav[j + i] -= tmp[i]; 
      outbv[j + i] -= inbv[j + i]; 
    }
  }
}

/* out += in*X^a */
// Assumes that a is even
void trlwe_mul_by_xai_addto_comp_vaes(TRLWE out, TRLWE in, int a){
  const int N = out->b->N;
  torus_polynomial_mul_by_xai_addto(out->b, in->b, a);
  const uint64_t id = in->a[0]->coeffs[0];
  __m512i * outav = (__m512i *) out->a[0]->coeffs;
  __m512i * outbv = (__m512i *) out->b->coeffs;
  __m512i * inbv = (__m512i *) in->b->coeffs;

  __m512i idv =  _mm512_maskz_set1_epi64 (0x55, id);
  __m512i cntv = {0, 0, 0, 1, 0, 2, 0, 3};
  
  const __m512i modN_mask = {-1, N/2 - 1, -1, N/2 - 1, -1, N/2 - 1, -1, N/2 - 1};
  cntv = _mm512_add_epi64 (cntv, idv);
  __m512i tmp[4], tmp2;

  a &= ((N<<1) - 1); // a % 2N
  if (!a) return;
  if (a < N) {
    __m512i av = _mm512_maskz_set1_epi64 (0xAA, (N-a)>>1);
    av = _mm512_add_epi64 (cntv, av);
    int i;
    for (i = 0; i < a/8 - 4; i+=4) {
      aes_next_16(tmp, &av);
      for (size_t j = 0; j < 4; j++){
        outav[i + j] -= tmp[j];
      }
    }
    for (; i < a/8; i++){
      aes_next_4(&tmp2, &av);
      outav[i] -= tmp2;
    }
    
    cntv = _mm512_and_epi64 (av, modN_mask);
    aes_next_4(&tmp[0], &av);
    av = cntv;
    aes_next_4(&tmp[1], &av);
    av = _mm512_and_epi64 (av, modN_mask);
    outav[a/8] += _mm512_mask_blend_epi64 (0xFFU<<(a&0x7), -tmp[0], tmp[1]);
    for (size_t i = 1; i < 4; i++){
      aes_next_4(&tmp[0], &av);
      outav[a/8 + i] += tmp[0];
    }

    for (size_t i = a/8 + 4; i < N/8; i+=4){
      aes_next_16(tmp, &av);
      for (size_t j = 0; j < 4 && (i + j) < N/8; j++){
        outav[i + j] += tmp[j];
        // printf("%lu %lu %lu\n", i, j, i + j);
      }
    }
  }else{
    a -= N;
    __m512i av = _mm512_maskz_set1_epi64 (0xAA, (N-a)>>1);
    av = _mm512_add_epi64 (cntv, av);
    for (int i = 0; i < a/8; i+=4) {
      aes_next_16(tmp, &av);
      for (size_t j = 0; j < 4; j++){
        outav[i + j] += tmp[j];
        outbv[i + j] += inbv[i + j];
      }
    }
    
    aes_next_4(&tmp[0], &av);
    av = cntv;
    aes_next_4(&tmp[1], &av);
    outav[a/8] += _mm512_mask_blend_epi64 (0xFF>>(a&0x7), +tmp[0], -tmp[1]);
    outbv[a/8] += _mm512_mask_blend_epi64 (0xFF>>(a&0x7), +inbv[a/8], -inbv[a/8]);
    for (size_t i = 1; i < 4; i++){
      aes_next_4(tmp, &av);
      outav[a/8 + i] -= tmp[0];
      outbv[a/8 + i] -= inbv[a/8 + i];
    }

    for (int i = a/8 + 4; i < N/8; i+=4){
      aes_next_16(tmp, &av);
      for (size_t j = 0; j < 4; j++){
        outav[i + j] -= tmp[j];
        outbv[i + j] -= inbv[i + j];
      }
    }
  }
}

void trgsw_mul_by_xai_addto_comp_vaes(TRGSW out, TRGSW in, int a){
  const int l = in->l, k = in->samples[0]->k;
  for (size_t i = 0; i < (k+1)*l; i++){
    trlwe_mul_by_xai_addto_comp_vaes(out->samples[i], in->samples[i], a);
  }
}