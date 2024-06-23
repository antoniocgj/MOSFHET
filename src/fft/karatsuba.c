#include "mosfhet.h"
#define MAX_DEPTH 1
__uint128_t * kara_tmp_out[MAX_DEPTH + 1], * kara_tmp_in1[MAX_DEPTH + 1], * kara_tmp_in2[MAX_DEPTH + 1];

__uint128_t int128_mul(uint64_t in1, uint64_t in2){
  __uint128_t output;
  #ifndef PORTABLE_BUILD
  long long unsigned int * out = (long long unsigned int *) &output;
  out[0] = _mulx_u64(in1, in2, &out[1]);
  #else
  output = ((__uint128_t) in1) * ((__uint128_t) in2);
  #endif
  return output;
}

void _debug_print128(char * msg, __uint128_t x){
  uint64_t * x64 = (uint64_t *) &x;
  printf("%s: 0x%016lx%016lx\n", msg, x64[1], x64[0]);
}


void poly_add_u128(__uint128_t * out, __uint128_t * p, __uint128_t * p2, int size){
  for (size_t i = 0; i < size; i++){
    out[i] = p[i] + p2[i];
  }
}

void poly_mul_int64to128(__uint128_t * out, uint64_t * p, uint64_t * p2, int size){
  const int new_size = size*2 - 1;
  memset(out, 0, sizeof(__uint128_t) * new_size);
  for (size_t i = 0; i < size; i++){
    for (size_t j = 0; j < size; j++){
      out[i+j] += int128_mul(p[i], p2[j]);
    }
  }
}

void poly_mul_u128(__uint128_t * out, __uint128_t * p, __uint128_t * p2, int size){
  const int new_size = size*2 - 1;
  memset(out, 0, sizeof(__uint128_t) * new_size);
  for (size_t i = 0; i < size; i++){
    for (size_t j = 0; j < size; j++){
      out[i+j] += p[i] * p2[j];
    }
  }
}

int init = 0;
void init_karatsuba(int size){ 
  if(init++) return;
  int m = size;
  for (size_t i = 0; i <= MAX_DEPTH; i++){
    kara_tmp_out[i] = (__uint128_t *)safe_aligned_malloc(sizeof(__uint128_t)*(m*2 - 1));
    kara_tmp_in1[i] = (__uint128_t *)safe_aligned_malloc(sizeof(__uint128_t)*m);
    kara_tmp_in2[i] = (__uint128_t *)safe_aligned_malloc(sizeof(__uint128_t)*m);
    m /= 2;
  }
}
 

void _karatsuba(__uint128_t * out, __uint128_t * in1, __uint128_t * in2, int size, int depth){
  if(depth > MAX_DEPTH){
    poly_mul_u128(out, in1, in2, size);
    return;
  }
  const int m = size/2;
  __uint128_t * a0 = in1;
  __uint128_t * a1 = &in1[m];
  __uint128_t * b0 = in2;
  __uint128_t * b1 = &in2[m];

  memset(out, 0, sizeof(__uint128_t)*(size*2 - 1));

  _karatsuba(kara_tmp_out[depth], a0, b0, m, depth + 1); // p0
  for (size_t i = 0; i < 2*m - 1; i++){
    out[i] += kara_tmp_out[depth][i];
    out[i + m] -= kara_tmp_out[depth][i];
  }
  _karatsuba(kara_tmp_out[depth], a1, b1, m, depth + 1); // p2
  for (size_t i = 0; i < 2*m - 1; i++){
    out[i + 1*m] -= kara_tmp_out[depth][i];
    out[i + 2*m] += kara_tmp_out[depth][i];
  }
  poly_add_u128(kara_tmp_in1[depth], a0, a1, m);
  poly_add_u128(kara_tmp_in2[depth], b0, b1, m);
  _karatsuba(kara_tmp_out[depth], kara_tmp_in1[depth], kara_tmp_in2[depth], m, depth + 1); // p1
  for (size_t i = 0; i < 2*m - 1; i++){
    out[i + 1*m] += kara_tmp_out[depth][i];
  }
}

void karatsuba_u128_scale64(uint64_t * out, uint64_t * in1, uint64_t * in2, int size, int bit_scale){
  init_karatsuba(size);
  for (size_t i = 0; i < size; i++){
    kara_tmp_in1[0][i] = (__uint128_t) in1[i];
    kara_tmp_in2[0][i] = (__uint128_t) in2[i];
  }
  _karatsuba(kara_tmp_out[0], kara_tmp_in1[0], kara_tmp_in2[0], size, 1);
  for (size_t i = 0; i < size; i++){
    out[i] = (uint64_t)(kara_tmp_out[0][i]>>bit_scale) - (uint64_t)(kara_tmp_out[0][size + i]>>bit_scale);
  }
}