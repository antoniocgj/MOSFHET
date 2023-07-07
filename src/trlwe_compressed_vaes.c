#include "mosfhet.h"
// k = 1
// aes setup
#define ID_SIZE 16

void aes_prng(uint8_t *output, uint64_t outlen, const uint8_t *input,  uint64_t inlen);
void aes_prgn_next_16(__m512i * out, __m512i * cnt);
void aes_prgn_next_4(__m512i * out, __m512i * cnt);

// size is the number of 64-bit words to generate
void aes_full_expand(__m512i * out, uint64_t * id, uint64_t size){
  aes_prng((uint8_t *) out, size*8, (uint8_t *) id, ID_SIZE);
}

void double_exp_normalize(__m512i * v, const uint64_t size){
  const __m512i clear_exp = _mm512_set1_epi64(0x800FFFFFFFFFFFFF);
  const __m512i exp = _mm512_set1_epi64(((uint64_t)(1023 + 64)) << 52);

  for (size_t i = 0; i < size; i++){
    v[i] = _mm512_ternarylogic_epi64 (v[i], clear_exp, exp, 0xEA);
  }
}

TRLWE trlwe_alloc_new_compressed_sample(int k, int N){
  TRLWE res;
  res = (TRLWE) safe_malloc(sizeof(*res));
  res->a = (TorusPolynomial *) safe_malloc(sizeof(TorusPolynomial));
  res->a[0] = polynomial_new_torus_polynomial(ID_SIZE/sizeof(Torus));
  res->b = polynomial_new_torus_polynomial(N);
  res->k = k;
  return res;
}

TRLWE_DFT trlwe_alloc_new_compressed_DFT_sample(int k, int N){
  TRLWE_DFT res;
  res = (TRLWE_DFT) safe_malloc(sizeof(*res));
  res->a = (DFT_Polynomial *) safe_malloc(sizeof(DFT_Polynomial));
  res->a[0] = polynomial_new_DFT_polynomial(ID_SIZE/sizeof(double));
  res->b = polynomial_new_DFT_polynomial(N);
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
  uint64_t * seed = out->a[0]->coeffs;

  // add error
  generate_torus_normal_random_array(out->b->coeffs, key->sigma, N);

  // internal product
  uint64_t seed_i[ID_SIZE/sizeof(uint64_t)];
  memcpy(seed_i, seed, ID_SIZE);
  for (size_t i = 0; i < key->k; i++){
    aes_full_expand((__m512i *)p_tmp->coeffs, seed_i, N);
    polynomial_mul_addto_torus(out->b, p_tmp, key->s[i]);
    seed_i[ID_SIZE/sizeof(uint64_t) - 1] += N/2;
  }

  if(m != NULL){
    for (size_t i = 0; i < N; i++){
      out->b->coeffs[i] += m->coeffs[i];
    }
  }
  free_polynomial(p_tmp);
}

void trlwe_compressed_DFT_sample(TRLWE_DFT out, TorusPolynomial m, TRLWE_Key key){
  const int N = key->s[0]->N, byte_size = ID_SIZE;
  TorusPolynomial p_tmp = polynomial_new_torus_polynomial(N);
  TorusPolynomial acc = polynomial_new_torus_polynomial(N);
  DFT_Polynomial p_tmp_dft = polynomial_new_DFT_polynomial(N);

  generate_random_bytes(byte_size, (uint8_t *) out->a[0]->coeffs);
  uint64_t * seed = (uint64_t *) out->a[0]->coeffs;

  // add error
  generate_torus_normal_random_array(acc->coeffs, key->sigma, N);

  // internal product
  uint64_t seed_i[ID_SIZE/sizeof(uint64_t)];
  memcpy(seed_i, seed, ID_SIZE);
  for (size_t i = 0; i < key->k; i++){
    aes_full_expand((__m512i *)p_tmp_dft->coeffs, seed_i, N);
    double_exp_normalize((__m512i *) p_tmp_dft->coeffs, N/8);
    polynomial_DFT_to_torus(p_tmp, p_tmp_dft);
    polynomial_mul_addto_torus(acc, p_tmp, key->s[i]);
    seed_i[ID_SIZE/sizeof(uint64_t) - 1] += N/2;
  }

  if(m != NULL){
    for (size_t i = 0; i < N; i++){
      acc->coeffs[i] += m->coeffs[i];
    }
  }

  polynomial_torus_to_DFT(out->b, acc);

  free_DFT_polynomial(p_tmp_dft);
  free_polynomial(p_tmp);
  free_polynomial(acc);
}

TRLWE trlwe_new_compressed_sample(TorusPolynomial m, TRLWE_Key key){
  const int N = key->s[0]->N;
  TRLWE res = trlwe_alloc_new_compressed_sample(key->k, N);
  trlwe_compressed_sample(res, m, key);
  return res;
}

TRLWE_DFT trlwe_new_compressed_DFT_sample(TorusPolynomial m, TRLWE_Key key){
  const int N = key->s[0]->N;
  TRLWE_DFT res = trlwe_alloc_new_compressed_DFT_sample(key->k, N);
  trlwe_compressed_DFT_sample(res, m, key);
  return res;
}


void trlwe_compressed_subto(TRLWE out, TRLWE in){
  const int N = out->b->N;

  const __m128i * id = (__m128i *) in->a[0]->coeffs;
  __m512i * outav = (__m512i *) out->a[0]->coeffs;
  __m512i * outbv = (__m512i *) out->b->coeffs;
  __m512i * inbv = (__m512i *) in->b->coeffs;

  __m128i cnt = _mm_loadu_si128(id);
  __m512i cntv = _mm512_broadcast_i64x2 (cnt);
  __m512i cntv2 = {0, 0, 0, 1, 0, 2, 0, 3};
  cntv = _mm512_add_epi64 (cntv, cntv2);
  __m512i tmp[4];

  for (size_t j = 0; j < N/8; j+=4){
    aes_prgn_next_16(tmp, &cntv);
    for (size_t i = 0; i < 4; i++){
      outav[j + i] -= tmp[i]; 
      outbv[j + i] -= inbv[j + i]; 
    }
  }
}

void trlwe_compressed_DFT_mul_addto(TRLWE_DFT out, DFT_Polynomial in1, TRLWE_DFT in2){
  const int N = out->b->N;

  // vector pointers
  const __m128i * id = (__m128i *) in2->a[0]->coeffs;
  __m512d * inbv = (__m512d *) in2->b->coeffs;
  __m512d * p = (__m512d *) in1->coeffs;
  __m512d * outav = (__m512d *) out->a[0]->coeffs;
  __m512d * outbv = (__m512d *) out->b->coeffs;

  // setup AES counters
  const __m128i cnt = _mm_loadu_si128(id);
  const __m512i cntv = _mm512_broadcast_i64x2 (cnt);
  const __m512i cntv2 = {0, 0, 0, 1, 0, N/4 + 0, 0, N/4 + 1};
  __m512i cntv_re = _mm512_add_epi64 (cntv, cntv2);
  __m512d tmp_re[4];
  __m512d * tmp_im = &tmp_re[2];

  // mul
  for (size_t j = 0; j < N/16; j+=2){
    aes_prgn_next_16((__m512i *)tmp_re, &cntv_re);
    double_exp_normalize((__m512i *)tmp_re, 4);
    for (size_t i = 0; i < 2; i++){
      const __m512d p_re = p[j + i];
      const __m512d p_im = p[j + i + N/16];
      const __m512d inbv_re = inbv[j + i];
      const __m512d inbv_im = inbv[j + i + N/16];

      // a
      const __m512d _1 = _mm512_fmsub_pd (p_im, tmp_im[i], outav[j + i]);
      const __m512d _2 = _mm512_fmadd_pd (p_im, tmp_re[i], outav[j + i + N/16]);
      outav[j + i] = _mm512_fmsub_pd (p_re, tmp_re[i],  _1);
      outav[j + i + N/16] = _mm512_fmadd_pd (p_re, tmp_im[i], _2);
      // b
      const __m512d _1b = _mm512_fmsub_pd (p_im, inbv_im, outbv[j + i]);
      const __m512d _2b = _mm512_fmadd_pd (p_im, inbv_re, outbv[j + i + N/16]);
      outbv[j + i] = _mm512_fmsub_pd (p_re, inbv_re,  _1b);
      outbv[j + i + N/16] = _mm512_fmadd_pd (p_re, inbv_im, _2b);
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
      aes_prgn_next_16(tmp, &av);
      for (size_t j = 0; j < 4; j++){
        outav[i + j] -= tmp[j];
      }
    }
    for (; i < a/8; i++){
      aes_prgn_next_4(&tmp2, &av);
      outav[i] -= tmp2;
    }
    
    cntv = _mm512_and_epi64 (av, modN_mask);
    aes_prgn_next_4(&tmp[0], &av);
    av = cntv;
    aes_prgn_next_4(&tmp[1], &av);
    av = _mm512_and_epi64 (av, modN_mask);
    outav[a/8] += _mm512_mask_blend_epi64 (0xFFU<<(a&0x7), -tmp[0], tmp[1]);
    for (size_t i = 1; i < 4; i++){
      aes_prgn_next_4(&tmp[0], &av);
      outav[a/8 + i] += tmp[0];
    }

    for (size_t i = a/8 + 4; i < N/8; i+=4){
      aes_prgn_next_16(tmp, &av);
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
      aes_prgn_next_16(tmp, &av);
      for (size_t j = 0; j < 4; j++){
        outav[i + j] += tmp[j];
        outbv[i + j] += inbv[i + j];
      }
    }
    
    aes_prgn_next_4(&tmp[0], &av);
    av = cntv;
    aes_prgn_next_4(&tmp[1], &av);
    outav[a/8] += _mm512_mask_blend_epi64 (0xFF>>(a&0x7), +tmp[0], -tmp[1]);
    outbv[a/8] += _mm512_mask_blend_epi64 (0xFF>>(a&0x7), +inbv[a/8], -inbv[a/8]);
    for (size_t i = 1; i < 4; i++){
      aes_prgn_next_4(tmp, &av);
      outav[a/8 + i] -= tmp[0];
      outbv[a/8 + i] -= inbv[a/8 + i];
    }

    for (int i = a/8 + 4; i < N/8; i+=4){
      aes_prgn_next_16(tmp, &av);
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