#include "mosfhet.h"

TRLWE trlwe_alloc_new_sample(int k, int N){
  TRLWE res;
  res = (TRLWE) safe_malloc(sizeof(*res));
  res->a = (TorusPolynomial *) safe_malloc(sizeof(TorusPolynomial) * k);
  for (size_t i = 0; i < k; i++){
    res->a[i] = polynomial_new_torus_polynomial(N);
  }
  res->b = polynomial_new_torus_polynomial(N);
  res->k = k;
  return res;
}

TRLWE * trlwe_alloc_new_sample_array(int count, int k, int N){
  TRLWE * res;
  res = (TRLWE *) safe_malloc(sizeof(TRLWE) * count);
  for (size_t i = 0; i < count; i++){
    res[i] = trlwe_alloc_new_sample(k, N);
  }
  return res;
}

void trlwe_save_sample(FILE * fd, TRLWE c){
  for (size_t i = 0; i < c->k; i++){
    fwrite(c->a[i]->coeffs, sizeof(Torus), c->a[0]->N, fd);
  }
  fwrite(c->b->coeffs, sizeof(Torus), c->b->N, fd);
}

void trlwe_load_sample(FILE * fd, TRLWE c){
  const int k = c->k, N = c->b->N;
  for (size_t i = 0; i < k; i++){
    fread(c->a[i]->coeffs, sizeof(Torus), N, fd);
  }
  fread(c->b->coeffs, sizeof(Torus), N, fd);
}

TRLWE trlwe_load_new_sample(FILE * fd, int k, int N){
  TRLWE res = trlwe_alloc_new_sample(k, N);
  trlwe_load_sample(fd, res);
  return res;
}

TRLWE_DFT * trlwe_alloc_new_DFT_sample_array(int count, int k, int N){
  TRLWE_DFT * res;
  res = (TRLWE_DFT *) safe_malloc(sizeof(TRLWE_DFT) * count);
  for (size_t i = 0; i < count; i++){
    res[i] = trlwe_alloc_new_DFT_sample(k, N);
  }
  return res;
}

TRLWE_DFT trlwe_alloc_new_DFT_sample(int k, int N){
  TRLWE_DFT res;
  res = (TRLWE_DFT) safe_malloc(sizeof(*res));
  res->a = (DFT_Polynomial *) safe_malloc(sizeof(DFT_Polynomial) * k);
  for (size_t i = 0; i < k; i++){
    res->a[i] = polynomial_new_DFT_polynomial(N);
  }
  res->b = polynomial_new_DFT_polynomial(N);
  res->k = k;
  return res;
}

void trlwe_save_DFT_sample(FILE * fd, TRLWE_DFT c){
  for (size_t i = 0; i < c->k; i++){
    fwrite(c->a[0]->coeffs, sizeof(double), c->a[0]->N, fd);
  }
  fwrite(c->b->coeffs, sizeof(double), c->b->N, fd);
}

TRLWE_DFT trlwe_load_new_DFT_sample(FILE * fd, int k, int N){
  TRLWE_DFT res = trlwe_alloc_new_DFT_sample(k, N);
  trlwe_load_DFT_sample(fd, res);
  return res;
}

void trlwe_load_DFT_sample(FILE * fd, TRLWE_DFT c){
  const int k = c->k, N = c->b->N;
  for (size_t i = 0; i < k; i++){
    fread(c->a[0]->coeffs, sizeof(double), N, fd);
  }
  fread(c->b->coeffs, sizeof(double), N, fd);
}

void free_trlwe(void * p_v){
  const TRLWE p = (TRLWE) p_v;
  for (size_t i = 0; i < p->k; i++){
    free_polynomial(p->a[i]);
  }
  free_polynomial(p->b);
  free(p->a);
  free(p);
}

void free_trlwe_array(void * p_v, int count){
  for (size_t i = 0; i < count; i++){
    free_trlwe(((void **)p_v)[i]);
  } 
  free(p_v);
}

TRLWE_Key trlwe_alloc_key(int N, int k, double sigma){
  TRLWE_Key res;
  res = (TRLWE_Key) safe_malloc(sizeof(*res));
  res->k = k;
  res->sigma = sigma;
  res->s = (IntPolynomial*) safe_malloc(k*sizeof(IntPolynomial));
  for (size_t i = 0; i < k; i++){
    res->s[i] = polynomial_new_torus_polynomial(N);  
  }
  return res;
}

// bound must be a power of 2
TRLWE_Key trlwe_new_bounded_key(int N, int k, uint64_t bound, double sigma){
  TRLWE_Key res = trlwe_alloc_key(N, k, sigma);
  for (size_t i = 0; i < k; i++){
    generate_random_bytes(N*sizeof(Torus), (uint8_t *) res->s[i]->coeffs);
    for (size_t j = 0; j < N; j++){
      res->s[i]->coeffs[j] &= (bound - 1);
      res->s[i]->coeffs[j] -= (bound >> 1) - 1;
    }    
  }
  return res;
}

TRLWE_Key trlwe_new_binary_key(int N, int k, double sigma){
  return trlwe_new_bounded_key(N, k, 2, sigma);
}

void trlwe_save_key(FILE * fd, TRLWE_Key key){
  fwrite(&key->k, sizeof(int), 1, fd);
  fwrite(&key->s[0]->N, sizeof(int), 1, fd);
  fwrite(&key->sigma, sizeof(double), 1, fd);
  for (size_t i = 0; i < key->k; i++){
    fwrite(key->s[i]->coeffs, sizeof(Torus), key->s[0]->N, fd);
  }
}

TRLWE_Key trlwe_load_new_key(FILE * fd){
  int k, N;
  double sigma;
  fread(&k, sizeof(int), 1, fd);
  fread(&N, sizeof(int), 1, fd);
  fread(&sigma, sizeof(double), 1, fd);
  TRLWE_Key key = trlwe_alloc_key(N, k, sigma);
  for (size_t i = 0; i < key->k; i++){
    fread(key->s[i]->coeffs, sizeof(Torus), key->s[0]->N, fd);
  }
  return key;
}

void free_trlwe_key(TRLWE_Key key){
  for (size_t i = 0; i < key->k; i++){
    free_polynomial(key->s[i]);
  }
  free(key->s);
  free(key);
}

void trlwe_noiseless_trivial_sample(TRLWE out, TorusPolynomial m){
  for (size_t i = 0; i < out->k; i++){
    memset(out->a[i]->coeffs, 0, sizeof(Torus)* out->b->N);
  }
  if(m != NULL) memcpy(out->b->coeffs, m->coeffs, sizeof(Torus) * out->b->N);
  else memset(out->b->coeffs, 0, sizeof(Torus) * out->b->N);
}

TRLWE trlwe_new_noiseless_trivial_sample(TorusPolynomial m, int k, int N){
  TRLWE res = trlwe_alloc_new_sample(k, N);
  trlwe_noiseless_trivial_sample(res, m);
  return res;
}

void trlwe_noiseless_trivial_DFT_sample(TRLWE_DFT out, DFT_Polynomial m){
  for (size_t i = 0; i < out->k; i++){
    memset(out->a[i]->coeffs, 0, sizeof(double)* out->b->N);
  }
  if(m != NULL) memcpy(out->b->coeffs, m->coeffs, sizeof(double) * out->b->N);
  else memset(out->b->coeffs, 0, sizeof(double) * out->b->N);
}

TRLWE_DFT trlwe_new_noiseless_trivial_DFT_sample(DFT_Polynomial m, int k, int N){
  TRLWE_DFT res = trlwe_alloc_new_DFT_sample(k, N);
  trlwe_noiseless_trivial_DFT_sample(res, m);
  return res;
}

void trlwe_sample(TRLWE out, TorusPolynomial m, TRLWE_Key key){
  const int N = key->s[0]->N, byte_size = sizeof(Torus) * N;

  for (size_t i = 0; i < key->k; i++){
    generate_random_bytes(byte_size, (uint8_t *) out->a[i]->coeffs);
  }

  // add error
  generate_torus_normal_random_array(out->b->coeffs, key->sigma, N);

  // internal product
  for (size_t i = 0; i < key->k; i++){
    polynomial_mul_addto_torus(out->b, out->a[i], key->s[i]);
  }

  if(m != NULL){
    for (size_t i = 0; i < m->N; i++){
      out->b->coeffs[i] += m->coeffs[i];
    }
  }
}

TRLWE trlwe_new_sample(TorusPolynomial m, TRLWE_Key key){
  TRLWE res = trlwe_alloc_new_sample(key->k, key->s[0]->N);
  trlwe_sample(res, m, key);
  return res;
}

void trlwe_phase(TorusPolynomial out, TRLWE in, TRLWE_Key key){
  const int N = key->s[0]->N, k = in->k, byte_size = sizeof(Torus) * N;
  memset(out->coeffs, 0, byte_size);
  for (size_t i = 0; i < k; i++){
    polynomial_mul_addto_torus(out, in->a[i], key->s[i]);
  }
  polynomial_sub_torus_polynomials(out, in->b, out);
}

void trlwe_add(TRLWE out, TRLWE in1, TRLWE in2){
  for (size_t i = 0; i < in1->k; i++){
    polynomial_add_torus_polynomials(out->a[i], in1->a[i], in2->a[i]);
  }
  polynomial_add_torus_polynomials(out->b, in1->b, in2->b);
}

void trlwe_copy(TRLWE out, TRLWE in){
  for (size_t i = 0; i < in->k; i++){
    polynomial_copy_torus_polynomial(out->a[i], in->a[i]);
  }
  polynomial_copy_torus_polynomial(out->b, in->b);
}


void trlwe_negate(TRLWE out, TRLWE in){
  for (size_t i = 0; i < in->k; i++){
    polynomial_negate_torus_polynomial(out->a[i], in->a[i]);
  }
  polynomial_negate_torus_polynomial(out->b, in->b);
}

void trlwe_DFT_copy(TRLWE_DFT out, TRLWE_DFT in){
  for (size_t i = 0; i < in->k; i++){
    polynomial_copy_DFT_polynomial(out->a[i], in->a[i]);
  }
  polynomial_copy_DFT_polynomial(out->b, in->b);
}

void trlwe_addto(TRLWE out, TRLWE in){
  trlwe_add(out, out, in);
}

void trlwe_DFT_add(TRLWE_DFT out, TRLWE_DFT in1,  TRLWE_DFT in2){
  for (size_t i = 0; i < in1->k; i++){
    polynomial_add_DFT_polynomials(out->a[i], in1->a[i], in2->a[i]);
  }
  polynomial_add_DFT_polynomials(out->b, in1->b, in2->b);
}

void trlwe_DFT_sub(TRLWE_DFT out, TRLWE_DFT in1,  TRLWE_DFT in2){
  for (size_t i = 0; i < in1->k; i++){
    polynomial_sub_DFT_polynomials(out->a[i], in1->a[i], in2->a[i]);
  }
  polynomial_sub_DFT_polynomials(out->b, in1->b, in2->b);
}

void trlwe_DFT_addto(TRLWE_DFT out, TRLWE_DFT in){
  trlwe_DFT_add(out, out, in);
}

#ifndef AVX512_OPT

void trlwe_sub(TRLWE out, TRLWE in1, TRLWE in2){
  for (size_t i = 0; i < in1->k; i++){
    polynomial_sub_torus_polynomials(out->a[i], in1->a[i], in2->a[i]);
  }
  polynomial_sub_torus_polynomials(out->b, in1->b, in2->b);
}

#else 
void trlwe_sub(TRLWE out, TRLWE in1, TRLWE in2){
  __m512i * a1 = (__m512i *) in1->a[0]->coeffs;
  __m512i * b1 = (__m512i *) in2->a[0]->coeffs;
  __m512i * c1 = (__m512i *) out->a[0]->coeffs;
  __m512i * a2 = (__m512i *) in1->b->coeffs;
  __m512i * b2 = (__m512i *) in2->b->coeffs;
  __m512i * c2 = (__m512i *) out->b->coeffs;
  for (size_t i = 0; i < in2->b->N/8; i++){
    c1[i] = _mm512_sub_epi64(a1[i], b1[i]);
    c2[i] = _mm512_sub_epi64(a2[i], b2[i]);
  }
}
#endif
void trlwe_subto(TRLWE out, TRLWE in){
  trlwe_sub(out, out, in);
}

void trlwe_DFT_mul_by_polynomial(TRLWE_DFT out, TRLWE_DFT in, DFT_Polynomial in2){
  const int k = in->k;
  polynomial_mul_DFT(out->a[0], in->a[0], in2);
  for (size_t i = 1; i < k; i++){
    polynomial_mul_addto_DFT(out->a[i], in->a[i], in2);
  }
  polynomial_mul_DFT(out->b, in->b, in2);
}

void trlwe_DFT_mul_addto_by_polynomial(TRLWE_DFT out, TRLWE_DFT in, DFT_Polynomial in2){
  const int k = in->k;
  polynomial_mul_addto_DFT(out->a[0], in->a[0], in2);
  for (size_t i = 1; i < k; i++){
    polynomial_mul_addto_DFT(out->a[i], in->a[i], in2);
  }
  polynomial_mul_addto_DFT(out->b, in->b, in2);
}

void trlwe_mul_by_xai(TRLWE out, TRLWE in, int a){
  const int k = in->k;
  for (size_t i = 0; i < k; i++){
    torus_polynomial_mul_by_xai(out->a[i], in->a[i], a);
  }
  torus_polynomial_mul_by_xai(out->b, in->b, a);
}

void trlwe_mul_by_xai_addto(TRLWE out, TRLWE in, int a){
  const int k = in->k;
  for (size_t i = 0; i < k; i++){
    torus_polynomial_mul_by_xai_addto(out->a[i], in->a[i], a);
  }
  torus_polynomial_mul_by_xai_addto(out->b, in->b, a);
}

void trlwe_mul_by_xai_minus_1(TRLWE out, TRLWE in, int a){
  const int k = in->k;
  for (size_t i = 0; i < k; i++){
    torus_polynomial_mul_by_xai_minus_1(out->a[i], in->a[i], a);
  }
  torus_polynomial_mul_by_xai_minus_1(out->b, in->b, a);
}

void trlwe_extract_tlwe_key(TLWE_Key out, TRLWE_Key in){
  const int N = in->s[0]->N, k = in->k;
  for (size_t i = 0; i < k; i++){
    for (size_t j = 0; j < N; j++){
      out->s[i*N + j] = in->s[i]->coeffs[j];
    }
  }
}

void trlwe_extract_tlwe(TLWE out, TRLWE in, int idx){
  const int N = in->b->N, k = in->k;
  for (size_t i = 0; i < k; i++){
    for (size_t j = 0; j <= idx; j++){
      out->a[i*N + j] = in->a[i]->coeffs[idx - j];
    }
    for (size_t j = idx + 1; j < N; j++){
      out->a[i*N + j] = -in->a[i]->coeffs[N + idx - j];
    }
  }
  out->b = in->b->coeffs[idx];
}

void trlwe_extract_tlwe_addto(TLWE out, TRLWE in, int idx){
  const int N = in->b->N, k = in->k;
  for (size_t i = 0; i < k; i++){
    for (size_t j = 0; j <= idx; j++){
      out->a[i*N + j] += in->a[i]->coeffs[idx - j];
    }
    for (size_t j = idx + 1; j < N; j++){
      out->a[i*N + j] += -in->a[i]->coeffs[N + idx - j];
    }
  }
  out->b += in->b->coeffs[idx];
}

void trlwe_extract_tlwe_subto(TLWE out, TRLWE in, int idx){
  const int N = in->b->N, k = in->k;
  for (size_t i = 0; i < k; i++){
    for (size_t j = 0; j <= idx; j++){
      out->a[i*N + j] -= in->a[i]->coeffs[idx - j];
    }
    for (size_t j = idx + 1; j < N; j++){
      out->a[i*N + j] -= -in->a[i]->coeffs[N + idx - j];
    }
  }
  out->b -= in->b->coeffs[idx];
}

void trlwe_mv_extract_tlwe(TLWE * out, TRLWE in, int amount){
  const int N = in->b->N;
  for (size_t i = 0; i < amount/2; i++){
    trlwe_extract_tlwe(out[i], in, i);
  }
  for (size_t i = amount/2; i < amount; i++){
    trlwe_extract_tlwe(out[i], in, N - 1 - (i - amount/2));
    tlwe_negate(out[i], out[i]);
  }
}

void trlwe_mv_extract_tlwe_scaling(TLWE out, TRLWE in, int scale){
  const int N = in->b->N, amount = scale;
  trlwe_extract_tlwe(out, in, amount/2);
  for (size_t i = amount/2 + 1; i < amount; i++){
    trlwe_extract_tlwe_subto(out, in, N - 1 - (i - amount/2));
  }
  for (size_t i = 0; i < amount/2; i++){
    trlwe_extract_tlwe_addto(out, in, i);
  }
}

void trlwe_mv_extract_tlwe_scaling_addto(TLWE out, TRLWE in, int scale){
  const int N = in->b->N, amount = scale;
  for (size_t i = amount/2; i < amount; i++){
    trlwe_extract_tlwe_subto(out, in, N - 1 - (i - amount/2));
  }
  for (size_t i = 0; i < amount/2; i++){
    trlwe_extract_tlwe_addto(out, in, i);
  }
}

void trlwe_mv_extract_tlwe_scaling_subto(TLWE out, TRLWE in, int scale){
  const int N = in->b->N, amount = scale;
  for (size_t i = amount/2; i < amount; i++){
    trlwe_extract_tlwe_addto(out, in, N - 1 - (i - amount/2));
  }
  for (size_t i = 0; i < amount/2; i++){
    trlwe_extract_tlwe_subto(out, in, i);
  }
}

void trlwe_to_DFT(TRLWE_DFT out, TRLWE in){
  for (size_t i = 0; i < in->k; i++){
    polynomial_torus_to_DFT(out->a[i], in->a[i]);
  }
  polynomial_torus_to_DFT(out->b, in->b);
}

void trlwe_from_DFT(TRLWE out, TRLWE_DFT in){
  for (size_t i = 0; i < in->k; i++){
    polynomial_DFT_to_torus(out->a[i], in->a[i]);
  }
  polynomial_DFT_to_torus(out->b, in->b);
}

void trlwe_decompose(TorusPolynomial * out, TRLWE in, int Bg_bit, int l){
  const int k = in->k, N = in->b->N;
  const uint64_t half_Bg = (1UL << (Bg_bit - 1));
  const uint64_t h_mask = (1UL << Bg_bit) - 1;
  const uint64_t word_size = sizeof(Torus)*8;

  uint64_t offset = 0;
  for (size_t i = 0; i < l; i++){
    offset += (1UL << (word_size - i * Bg_bit - 1));
  }
  
  for (size_t i = 0; i < l; i++) {
    const uint64_t h_bit = word_size - (i + 1) * Bg_bit;
    for (size_t j = 0; j < k; j++){
      for (size_t c = 0; c < N; c++){
        const uint64_t coeff_off = in->a[j]->coeffs[c] + offset;
        out[j*l + i]->coeffs[c] = ((coeff_off>>h_bit) & h_mask) - half_Bg;
      }
    }
    for (size_t c = 0; c < N; c++){
      const uint64_t coeff_off = in->b->coeffs[c] + offset;
      out[k*l + i]->coeffs[c] = ((coeff_off>>h_bit) & h_mask) - half_Bg;
    }
  }
}

void trlwe_torus_packing(TRLWE out, Torus * in, int size){
  trlwe_noiseless_trivial_sample(out, 0);
  for (size_t i = 0; i < out->b->N; i++){
    out->b->coeffs[i] = in[i/(out->b->N/size)];
  }
}

/* */
void trlwe_torus_packing_many_LUT(TRLWE out, Torus * in, int lut_size, int n_luts){
  trlwe_noiseless_trivial_sample(out, 0);
  for (size_t i = 0; i < lut_size; i++){
    for (size_t j = 0; j < n_luts; j++){
      for (size_t k = 0; k < out->b->N/(lut_size*n_luts); k++){
        out->b->coeffs[(i*n_luts + j)*out->b->N/(lut_size*n_luts) + k] = in[j*lut_size + i];
      }
    }
  }
}



// TRLWE tensor product (WIP) 
void trlwe_tensor_prod(TRLWE out, TRLWE in1, TRLWE in2, int precision, TRLWE_KS_Key rl_key){
  const int N = in1->b->N, bit_size = sizeof(Torus)*8;
  assert(in1->k == 1 && in2->k == 1);
  TorusPolynomial tmp = polynomial_new_torus_polynomial(N);
  TRLWE t = trlwe_alloc_new_sample(1, N);
  // T = A1 * A2
  polynomial_full_mul_with_scale(t->a[0], in1->a[0], in2->a[0], bit_size, bit_size - precision);
  for (size_t i = 0; i < N; i++) t->b->coeffs[i] = 0;
  // A = A1*B2 + B1*A2
  polynomial_full_mul_with_scale(out->a[0], in1->a[0], in2->b, bit_size, bit_size - precision);
  polynomial_full_mul_with_scale(tmp, in1->b, in2->a[0], bit_size, bit_size - precision);
  polynomial_addto_torus_polynomial(out->a[0], tmp);
  // B = B1*B2
  polynomial_full_mul_with_scale(out->b, in1->b, in2->b, bit_size, bit_size - precision);
  // Relinearization
  trlwe_keyswitch(t, t, rl_key);
  trlwe_subto(out, t);
  // Free
  free_polynomial(tmp);
  free_trlwe(t);
}

void trlwe_tensor_prod_FFT(TRLWE out, TRLWE in1, TRLWE in2, int precision, TRLWE_KS_Key rl_key){
  const int N = in1->b->N, half_prec = sizeof(Torus)*8 - (sizeof(Torus)*8 - precision)/2;
  assert(in1->k == 1 && in2->k == 1);
  TorusPolynomial tmp = polynomial_new_torus_polynomial(N);
  DFT_Polynomial tmp_DFT = polynomial_new_DFT_polynomial(N);
  TRLWE_DFT t = trlwe_alloc_new_DFT_sample(1, N);
  TRLWE t2 = trlwe_alloc_new_sample(1, N);
  DFT_Polynomial A1 = polynomial_new_DFT_polynomial(N);
  DFT_Polynomial A2 = polynomial_new_DFT_polynomial(N);
  DFT_Polynomial B1 = polynomial_new_DFT_polynomial(N);
  DFT_Polynomial B2 = polynomial_new_DFT_polynomial(N);
  // T = A1 * A2
  polynomial_torus_scale(tmp, in1->a[0], half_prec);
  polynomial_torus_to_DFT(A1, tmp);
  polynomial_torus_scale(tmp, in2->a[0], half_prec);
  polynomial_torus_to_DFT(A2, tmp);
  polynomial_mul_DFT(t->a[0], A1, A2);
  for (size_t i = 0; i < N; i++) t->b->coeffs[i] = 0.;
  // A = A1*B2 + B1*A2
  polynomial_torus_scale(tmp, in1->b, half_prec);
  polynomial_torus_to_DFT(B1, tmp);
  polynomial_torus_scale(tmp, in2->b, half_prec);
  polynomial_torus_to_DFT(B2, tmp);
  polynomial_mul_DFT(tmp_DFT, A1, B2);
  polynomial_mul_addto_DFT(tmp_DFT, B1, A2);
  polynomial_DFT_to_torus(out->a[0], tmp_DFT);
  // B = B1*B2
  polynomial_mul_DFT(tmp_DFT, B1, B2);
  polynomial_DFT_to_torus(out->b, tmp_DFT);
  // Relinearization
  trlwe_from_DFT(t2, t);
  trlwe_keyswitch(t2, t2, rl_key);
  trlwe_subto(out, t2);
  // Free
  free_polynomial(tmp);
  free_polynomial(tmp_DFT);
  free_trlwe(t);
  free_trlwe(t2);
  free_polynomial(A1);
  free_polynomial(A2);
  free_polynomial(B1);
  free_polynomial(B2);
}


// EvalAuto
void trlwe_eval_automorphism(TRLWE out, TRLWE in, uint64_t gen, TRLWE_KS_Key ks_key){
  polynomial_permute(out->a[0], in->a[0], gen);
  polynomial_permute(out->b, in->b, gen);
  trlwe_keyswitch(out, out, ks_key);
}
