#include "mosfhet.h"

TorusPolynomial polynomial_new_torus_polynomial(int N){
  TorusPolynomial res;
  res = (TorusPolynomial) safe_malloc(sizeof(*res));
  res->coeffs = (Torus *) safe_aligned_malloc(sizeof(Torus) * N);
  res->N = N;
  return res;
}

BinaryPolynomial polynomial_new_binary_polynomial(int N){
  BinaryPolynomial res;
  res = (BinaryPolynomial) safe_malloc(sizeof(*res));
  res->coeffs = (Binary *) safe_malloc(sizeof(Binary) * N);
  res->N = N;
  return res;
}

TorusPolynomial * polynomial_new_array_of_torus_polynomials(int N, int size){
  TorusPolynomial * res = (TorusPolynomial *) safe_malloc(sizeof(TorusPolynomial)*size);
  for (size_t i = 0; i < size; i++) res[i] = polynomial_new_torus_polynomial(N);
  return res;
}


DFT_Polynomial polynomial_new_DFT_polynomial(int N){
  DFT_Polynomial res;
  res = (DFT_Polynomial) safe_malloc(sizeof(*res));
  res->coeffs = (double *) safe_aligned_malloc(sizeof(double) * N);
  res->N = N;
  return res;
}

DFT_Polynomial * polynomial_new_array_of_polynomials_DFT(int N, int size){
  DFT_Polynomial * res = (DFT_Polynomial *) safe_malloc(sizeof(DFT_Polynomial)*size);
  for (size_t i = 0; i < size; i++) res[i] = polynomial_new_DFT_polynomial(N);
  return res;
}

void free_polynomial(void * p){
  free(((TorusPolynomial) p)->coeffs);
  free(p);
}

void free_DFT_polynomial(DFT_Polynomial p){
  free(p->coeffs);
  free(p);
}

void free_array_of_polynomials(void * p, int size){
  for (size_t i = 0; i < size; i++) free_polynomial(((void **) p)[i]);
  free(p);
}

void polynomial_decompose(TorusPolynomial * out, TorusPolynomial in, int Bg_bit, int l){
  const int N = in->N, bit_size = sizeof(Torus)*8;
  const uint64_t half_Bg = (1UL << (Bg_bit - 1));
  const uint64_t h_mask = (1UL << Bg_bit) - 1;

  uint64_t offset = 0;
  for (size_t i = 0; i < l; i++){
    offset += (1UL << (bit_size - i * Bg_bit - 1));
  }
  
  for (size_t i = 0; i < l; i++) {
    const uint64_t h_bit = bit_size - (i + 1) * Bg_bit;
    for (size_t c = 0; c < N; c++){
      const uint64_t coeff_off = in->coeffs[c] + offset;
      out[i]->coeffs[c] = ((coeff_off>>h_bit) & h_mask) - half_Bg;
    }
  }
}

void polynomial_decompose_i(TorusPolynomial out, TorusPolynomial in, int Bg_bit, int l, int i){
  const int N = in->N, bit_size = sizeof(Torus)*8;
  const uint64_t half_Bg = (1ULL << (Bg_bit - 1));
  const uint64_t h_mask = (1ULL << Bg_bit) - 1;
  const uint64_t h_bit = bit_size - (i + 1) * Bg_bit;

  uint64_t offset = 1ULL << (bit_size - l * Bg_bit - 1);
  for (size_t i = 0; i < l; i++){
    offset += (1ULL << (bit_size - i * Bg_bit - 1));
  }
  
  for (size_t c = 0; c < N; c++){
    const uint64_t coeff_off = in->coeffs[c] + offset;
    out->coeffs[c] = ((coeff_off>>h_bit) & h_mask) - half_Bg;
  }
}

/* Polynomial Arithmetic */

/* out = in1 + in2 */
void polynomial_add_torus_polynomials(TorusPolynomial out, TorusPolynomial in1, TorusPolynomial in2){
  for (size_t i = 0; i < in2->N; i++){
    out->coeffs[i] = in1->coeffs[i] + in2->coeffs[i];
  }
}

void polynomial_add_DFT_polynomials(DFT_Polynomial out, DFT_Polynomial in1, DFT_Polynomial in2){
  for (size_t i = 0; i < in2->N; i++){
    out->coeffs[i] = in1->coeffs[i] + in2->coeffs[i];
  }
}

void polynomial_sub_DFT_polynomials(DFT_Polynomial out, DFT_Polynomial in1, DFT_Polynomial in2){
  for (size_t i = 0; i < in2->N; i++){
    out->coeffs[i] = in1->coeffs[i] - in2->coeffs[i];
  }
}

/* out = in */
void polynomial_copy_torus_polynomial(TorusPolynomial out, TorusPolynomial in){
  memcpy(out->coeffs, in->coeffs, in->N*sizeof(Torus));
}

/* out = -in */
void polynomial_negate_torus_polynomial(TorusPolynomial out, TorusPolynomial in){
  for (size_t i = 0; i < in->N; i++){
    out->coeffs[i] = -in->coeffs[i];
  }
}

/* out = in */
void polynomial_copy_DFT_polynomial(DFT_Polynomial out, DFT_Polynomial in){
  memcpy(out->coeffs, in->coeffs, in->N*sizeof(double));
}

/* out += in */
void polynomial_addto_torus_polynomial(TorusPolynomial out, TorusPolynomial in){
  for (size_t i = 0; i < in->N; i++){
    out->coeffs[i] += in->coeffs[i];
  }
}

/* out = in1 - in2 */
void polynomial_sub_torus_polynomials(TorusPolynomial out, TorusPolynomial in1, TorusPolynomial in2){
  #ifdef AVX512_OPT
  __m512i * a = (__m512i *) in1->coeffs;
  __m512i * b = (__m512i *) in2->coeffs;
  __m512i * c = (__m512i *) out->coeffs;
  for (size_t i = 0; i < in2->N/8; i++){
    c[i] = _mm512_sub_epi64(a[i], b[i]);
  }
  #else
  for (size_t i = 0; i < in2->N; i++){
    out->coeffs[i] = in1->coeffs[i] - in2->coeffs[i];
  }
  #endif
}

/* out -= in */
void polynomial_subto_torus_polynomial(TorusPolynomial out, TorusPolynomial in){
  polynomial_sub_torus_polynomials(out, out, in);
}

/* out = in*X^a */
void torus_polynomial_mul_by_xai(TorusPolynomial out, TorusPolynomial in, int a){
  const int N = out->N;
  a &= ((N<<1) - 1); // a % 2N
  if (!a) return;
  if (a < N) {
    for (int i = 0; i < a; i++) out->coeffs[i] = -in->coeffs[i - a + N];
    for (int i = a; i < N; i++) out->coeffs[i] = in->coeffs[i - a];
  }else{
    for (int i = 0; i < a - N; i++) out->coeffs[i] = in->coeffs[i - a + 2*N];
    for (int i = a - N; i < N; i++) out->coeffs[i] = -in->coeffs[i - a + N];
  }
}

/* out += in*X^a */
void torus_polynomial_mul_by_xai_addto(TorusPolynomial out, TorusPolynomial in, int a){
  const int N = out->N;
  a &= ((N<<1) - 1); // a % 2N
  if (!a) return;
  if (a < N) {
    for (int i = 0; i < a; i++) out->coeffs[i] += -in->coeffs[i - a + N];
    for (int i = a; i < N; i++) out->coeffs[i] += in->coeffs[i - a];
  }else{
    for (int i = 0; i < a - N; i++) out->coeffs[i] += in->coeffs[i - a + 2*N];
    for (int i = a - N; i < N; i++) out->coeffs[i] += -in->coeffs[i - a + N];
  }
}

/* out = in*(X^a - 1)*/
void torus_polynomial_mul_by_xai_minus_1(TorusPolynomial out, TorusPolynomial in, int a){
  const int N = out->N;
  a &= ((N<<1) - 1); // a % 2N
  if (!a) return;
  if (a < N) {
    for (int i = 0; i < a; i++) out->coeffs[i] = -in->coeffs[i - a + N] - in->coeffs[i];
    for (int i = a; i < N; i++) out->coeffs[i] = in->coeffs[i - a] - in->coeffs[i];
  }else{
    for (int i = 0; i < a - N; i++) out->coeffs[i] = in->coeffs[i - a + 2*N] - in->coeffs[i];
    for (int i = a - N; i < N; i++) out->coeffs[i] = -in->coeffs[i - a + N] - in->coeffs[i];
  }
}

/* out += in1*in2 */
void polynomial_naive_mul_addto_torus_binary(TorusPolynomial out, TorusPolynomial in1, BinaryPolynomial in2){
  const int N = in2->N;
  for (size_t i = 0; i < N; i++){
    if(in2->coeffs[i]){
      for (size_t j = i; j < N; j++){
        out->coeffs[j] += in1->coeffs[j - i];
      }
      for (size_t j = 0; j < i; j++){
        out->coeffs[j] -= in1->coeffs[N + j - i];
      }
    }
  }
}

/* out += in1*in2 */
void polynomial_naive_mul_addto_torus(TorusPolynomial out, TorusPolynomial in1, TorusPolynomial in2){
  const int N = in2->N;
  for (size_t i = 0; i < N; i++){
    for (size_t j = i; j < N; j++){
      out->coeffs[j] += in1->coeffs[j - i]*in2->coeffs[i];
    }
    for (size_t j = 0; j < i; j++){
      out->coeffs[j] -= in1->coeffs[N + j - i]*in2->coeffs[i];
    }
  }
}


void polynomial_mul_torus(TorusPolynomial out, TorusPolynomial in1, TorusPolynomial in2){
  const int N = in2->N;
  // alloc temporaries
  static DFT_Polynomial * tmp_pool[32] = {NULL}; 
  if(tmp_pool[N>>10] == NULL) tmp_pool[N>>10] = polynomial_new_array_of_polynomials_DFT(N, 3);
  DFT_Polynomial * tmp = tmp_pool[N>>10];
  // dft mul
  polynomial_torus_to_DFT(tmp[1], in1);
  polynomial_torus_to_DFT(tmp[2], in2);
  polynomial_mul_DFT(tmp[0], tmp[1], tmp[2]);
  polynomial_DFT_to_torus(out, tmp[0]);
}

void polynomial_mul_addto_torus(TorusPolynomial out, TorusPolynomial in1, TorusPolynomial in2){
  const int N = in2->N;
  // alloc temporaries
  static TorusPolynomial tmp_pool[32] = {NULL}; 
  if(tmp_pool[N>>10] == NULL) tmp_pool[N>>10] = polynomial_new_torus_polynomial(N);
  TorusPolynomial tmp = tmp_pool[N>>10];
  polynomial_mul_torus(tmp, in1, in2);
  polynomial_addto_torus_polynomial(out, tmp);
}

/* out = in1*in2 */
void polynomial_naive_mul_torus(TorusPolynomial out, TorusPolynomial in1, TorusPolynomial in2){
  const int N = in2->N;
  for (size_t j = 0; j < N; j++){
    out->coeffs[j] = in1->coeffs[j]*in2->coeffs[0];
  }
  for (size_t i = 1; i < N; i++){
    for (size_t j = i; j < N; j++){
      out->coeffs[j] += in1->coeffs[j - i]*in2->coeffs[i];
    }
    for (size_t j = 0; j < i; j++){
      out->coeffs[j] -= in1->coeffs[N + j - i]*in2->coeffs[i];
    }
  }
}

/* out = in1*in2 */
void polynomial_naive_mul_binary(BinaryPolynomial out, BinaryPolynomial in1, BinaryPolynomial in2){
  const int N = in2->N;
  for (size_t j = 0; j < N; j++){
    out->coeffs[j] = in1->coeffs[j]*in2->coeffs[0];
  }
  for (size_t i = 1; i < N; i++){
    for (size_t j = i; j < N; j++){
      out->coeffs[j] += in1->coeffs[j - i]*in2->coeffs[i];
    }
    for (size_t j = 0; j < i; j++){
      out->coeffs[j] -= in1->coeffs[N + j - i]*in2->coeffs[i];
    }
  }
}

/* out[i] = round(in[i] * 2^log_scale), for i \\in [0,N) */
void polynomial_torus_scale(TorusPolynomial out, TorusPolynomial in, int log_scale){
  for (size_t i = 0; i < in->N; i++){
    out->coeffs[i] = (Torus) torus2int(in->coeffs[i], log_scale);
  }
}


#if defined(USE_SPQLIOS)
#include "../src/fft/spqlios/spqlios-fft.h"
FFT_Processor_Spqlios fft_proc[32] = {NULL};

void init_fft(int N){
  if(!fft_proc[N >> 10]) fft_proc[N >> 10] = new_FFT_Processor_Spqlios(N);
}
#else
#include "./fft/ffnt/ffnt.h"
FFT_Processor_FFNT fft_proc[8] = {NULL};

void init_fft(int N){
  if(!fft_proc[N >> 10]) fft_proc[N >> 10] = new_FFT_Processor_FFNT(N);
}
#endif

void polynomial_DFT_to_torus(TorusPolynomial out, const DFT_Polynomial in){
#ifdef TORUS32
  execute_direct_torus32(out->coeffs, in->coeffs, fft_proc[in->N >> 10]);
#else
  execute_direct_torus64(out->coeffs, in->coeffs, fft_proc[in->N >> 10]);
#endif
}

void polynomial_torus_to_DFT(DFT_Polynomial out, TorusPolynomial in){
  init_fft(in->N);
#ifdef TORUS32
  execute_reverse_torus32(out->coeffs, in->coeffs, fft_proc[in->N >> 10]);
#else
  execute_reverse_torus64(out->coeffs, in->coeffs, fft_proc[in->N >> 10]);
#endif
}


/* out = in1*in2 */
void polynomial_mul_DFT(DFT_Polynomial out, DFT_Polynomial in1, DFT_Polynomial in2){
  const int N = in1->N;
  #ifdef AVX512_OPT
  __m512d * a = (__m512d *) in1->coeffs;
  __m512d * b = (__m512d *) in2->coeffs;
  __m512d * c = (__m512d *) out->coeffs;
  for (int i = 0; i < N / 16; i++) {
    const __m512d a_im = a[i + N/16];
    const __m512d a_re = a[i];
    const __m512d b_re = b[i];
    const __m512d b_im = b[i + N/16];
    const __m512d _1 = _mm512_mul_pd (a_im, b_im);
    c[i] = _mm512_fmsub_pd (a_re, b_re,  _1);
    const __m512d _2 = _mm512_mul_pd (a_im, b_re);
    c[i + N/16] = _mm512_fmadd_pd (a_re, b_im,  _2);
  }
  #else
  for (int i = 0; i < N / 2; i++) {
    out->coeffs[i] = in1->coeffs[i] * in2->coeffs[i] - in1->coeffs[i + N / 2] * in2->coeffs[i + N / 2];
    out->coeffs[i + N / 2] = in1->coeffs[i + N / 2] * in2->coeffs[i] + in1->coeffs[i] * in2->coeffs[i + N / 2];
  }
  #endif
}



/* out += in1*in2 */
void polynomial_mul_addto_DFT(DFT_Polynomial out, DFT_Polynomial in1, DFT_Polynomial in2){
  const int N = in1->N;
  #ifdef AVX512_OPT
  __m512d * a = (__m512d *) in1->coeffs;
  __m512d * b = (__m512d *) in2->coeffs;
  __m512d * c = (__m512d *) out->coeffs;
  for (int i = 0; i < N / 16; i++) {
    const __m512d _1 = _mm512_fmsub_pd (a[i + N/16], b[i + N/16], c[i]);
    c[i] = _mm512_fmsub_pd (a[i], b[i],  _1);
    const __m512d _2 = _mm512_fmadd_pd (a[i + N/16], b[i], c[i + N/16]);
    c[i + N/16] = _mm512_fmadd_pd (a[i], b[i + N/16],  _2);
  }
  #else
  for (int i = 0; i < N / 2; i++) {
    out->coeffs[i] += in1->coeffs[i] * in2->coeffs[i] - in1->coeffs[i + N / 2] * in2->coeffs[i + N / 2];
    out->coeffs[i + N / 2] += in1->coeffs[i + N / 2] * in2->coeffs[i] + in1->coeffs[i] * in2->coeffs[i + N / 2];
  }
  #endif
}

void karatsuba_u128_scale64(uint64_t * out, uint64_t * in1, uint64_t * in2, int size, int bit_scale);
void polynomial_full_mul_with_scale(TorusPolynomial out, TorusPolynomial in1, TorusPolynomial in2, int bit_size, int bit_scale){
#ifndef TORUS32
  const int N = in1->N;
  assert(N == in2->N && N == out->N);
  karatsuba_u128_scale64(out->coeffs, in1->coeffs, in2->coeffs, N, bit_scale);
#else
  // TODO: use FFT (do we have a use for this function with TORUS32?)
#endif
}

// Galois transform: x^i to x^(gen*i)
// Code adapted from Lattigo 
// Src: https://github.com/tuneinsight/lattigo/blob/master/ring/ring_automorphism.go, Apache 2 License
void polynomial_permute(TorusPolynomial restrict out, TorusPolynomial restrict in, uint64_t gen){
  const uint64_t N = in->N, mask = N - 1;
  assert(out != in);
  for (size_t i = 0; i < N; i++){
    const uint64_t idx = i*gen;
    if(idx&N) out->coeffs[idx&mask] = -in->coeffs[i];
    else out->coeffs[idx&mask] = in->coeffs[i];
  }
}