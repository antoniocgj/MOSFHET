#include "mosfhet.h"

void free_trgsw(void * p_v){
  TRGSW p = (TRGSW) p_v;
  const int l = p->l, k  = p->samples[0]->k;
  for (size_t i = 0; i < l * (k + 1); i++){
    free_trlwe(p->samples[i]);
  }
  free(p->samples);
  free(p);
}

void free_trgsw_array(void * p_v, int count){
  for (size_t i = 0; i < count; i++){
    free_trgsw(((void **) p_v)[i]);
  }
  free(p_v);
}

TRGSW_Key trgsw_new_key(TRLWE_Key trlwe_key, int l, int Bg_bit){
  TRGSW_Key res;
  res = (TRGSW_Key) safe_malloc(sizeof(*res));
  res->trlwe_key = trlwe_key;
  res->l = l;
  res->Bg_bit = Bg_bit;
  return res;
}

TRGSW_Key trgsw_load_new_key(FILE * fd){
  TRGSW_Key res;
  res = (TRGSW_Key) safe_malloc(sizeof(*res));
  fread(&res->l, sizeof(int), 1, fd);
  fread(&res->Bg_bit, sizeof(int), 1, fd);
  res->trlwe_key = trlwe_load_new_key(fd);
  return res;
}

void trgsw_save_key(FILE * fd, TRGSW_Key key){
  fwrite(&key->l, sizeof(int), 1, fd);
  fwrite(&key->Bg_bit, sizeof(int), 1, fd);
  trlwe_save_key(fd, key->trlwe_key);
}

void free_trgsw_key(TRGSW_Key key){
  free(key);
}

TRGSW trgsw_alloc_new_sample(int l, int Bg_bit, int k, int N){
  TRGSW res;
  res = (TRGSW) safe_malloc(sizeof(*res));
  res->samples = (TRLWE *) safe_malloc(sizeof(TRLWE) * l * (k + 1));
  for (size_t i = 0; i < l * (k + 1); i++){
    res->samples[i] = trlwe_alloc_new_sample(k, N);
  }
  res->Bg_bit = Bg_bit;
  res->l = l;
  return res;
}

void trgsw_save_sample(FILE * fd, TRGSW c){
  for (size_t i = 0; i < c->l * (c->samples[0]->k + 1); i++){
    trlwe_save_sample(fd, c->samples[i]);
  }
}

TRGSW trgsw_load_new_sample(FILE * fd, int l, int Bg_bit, int k, int N){
  TRGSW res = trgsw_alloc_new_sample(l, Bg_bit, k, N);
  for (size_t i = 0; i < l * (k + 1); i++){
    trlwe_load_sample(fd, res->samples[i]);
  }
  return res;
}

void trgsw_load_sample(FILE * fd, TRGSW c){
  for (size_t i = 0; i < c->l * (c->samples[0]->k + 1); i++){
    trlwe_load_sample(fd, c->samples[i]);
  }
}

void trgsw_save_DFT_sample(FILE * fd, TRGSW_DFT c){
  for (size_t i = 0; i < c->l * (c->samples[0]->k + 1); i++){
    trlwe_save_DFT_sample(fd, c->samples[i]);
  }
}

TRGSW_DFT trgsw_load_new_DFT_sample(FILE * fd, int l, int Bg_bit, int k, int N){
  TRGSW_DFT res = trgsw_alloc_new_DFT_sample(l, Bg_bit, k, N);
  for (size_t i = 0; i < l * (k + 1); i++){
    trlwe_load_DFT_sample(fd, res->samples[i]);
  }
  return res;
}

void trgsw_load_DFT_sample(FILE * fd, TRGSW_DFT out){
  for (size_t i = 0; i < out->l * (out->samples[0]->k + 1); i++){
    trlwe_load_DFT_sample(fd, out->samples[i]);
  }
}

TRGSW * trgsw_alloc_new_sample_array(int count, int l, int Bg_bit, int k, int N){
  TRGSW * res;
  res = (TRGSW *) safe_malloc(sizeof(TRGSW)*count);
  for (size_t i = 0; i < count; i++){
    res[i] = trgsw_alloc_new_sample(l, Bg_bit, k, N);
  }
  return res;
}

TRGSW_DFT trgsw_alloc_new_DFT_sample(int l, int Bg_bit, int k, int N){
  TRGSW_DFT res;
  res = (TRGSW_DFT) safe_malloc(sizeof(*res));
  res->samples = (TRLWE_DFT *) safe_malloc(sizeof(TRLWE_DFT) * l * (k + 1));
  for (size_t i = 0; i < l * (k + 1); i++){
    res->samples[i] = trlwe_alloc_new_DFT_sample(k, N);
  }
  res->Bg_bit = Bg_bit;
  res->l = l;
  return res;
}

TRGSW_DFT * trgsw_alloc_new_DFT_sample_array(int count, int l, int Bg_bit, int k, int N){
  TRGSW_DFT * res;
  res = (TRGSW_DFT *) safe_malloc(sizeof(TRGSW_DFT)*count);
  for (size_t i = 0; i < count; i++){
    res[i] = trgsw_alloc_new_DFT_sample(l, Bg_bit, k, N);
  }
  return res;
}

void trgsw_noiseless_trivial_sample(TRGSW out, Torus m, int l, int Bg_bit, int k, int N){
  for (size_t i = 0; i < (k + 1)*l; i++){
    trlwe_noiseless_trivial_sample(out->samples[i], NULL);
  }
  
  for (size_t i = 0; i < l; i++) {
    const uint64_t h = 1UL << (sizeof(Torus)*8 - (i + 1) * Bg_bit);
    for (size_t j = 0; j < k; j++){
      out->samples[j*l + i]->a[j]->coeffs[0] += m * h;
    }
    out->samples[k*l + i]->b->coeffs[0] += m * h;
  }
}


TRGSW trgsw_new_noiseless_trivial_sample(Torus m, int l, int Bg_bit, int k, int N){
  TRGSW res = trgsw_alloc_new_sample(l, Bg_bit, k, N);
  trgsw_noiseless_trivial_sample(res, m, l, Bg_bit, k, N);
  return res;
}

/* TRGSW_key(mX^e) */
void trgsw_monomial_sample(TRGSW out, int64_t m, int e, TRGSW_Key key){
  const int l = key->l, k = key->trlwe_key->k, Bg_bit = key->Bg_bit, N = key->trlwe_key->s[0]->N;
  assert(out->l == key->l);
  if(e&N) m *= -1;
  e &= (N - 1);
  for (size_t i = 0; i < l * (k + 1); i++){
    trlwe_sample(out->samples[i], NULL, key->trlwe_key);
  }

  for (size_t i = 0; i < l; i++) {
    const Torus h = 1UL << (sizeof(Torus)*8 - (i + 1) * Bg_bit);
    for (size_t j = 0; j < k; j++){
      out->samples[j*l + i]->a[j]->coeffs[e] += m * h;
    }
    out->samples[k*l + i]->b->coeffs[e] += m * h;
  }
}

void trgsw_monomial_DFT_sample(TRGSW_DFT out, int64_t m, int e, TRGSW_Key key){
  TRGSW tmp = trgsw_alloc_new_sample(out->l, out->Bg_bit, out->samples[0]->k, out->samples[0]->b->N);
  trgsw_monomial_sample(tmp, m, e, key);
  trgsw_to_DFT(out, tmp);
  free_trgsw(tmp);
}

/* TRGSW_key(mX^e) */
TRGSW trgsw_new_monomial_sample(int64_t m, int e, TRGSW_Key key){
  const int l = key->l, k = key->trlwe_key->k, Bg_bit = key->Bg_bit, N = key->trlwe_key->s[0]->N;
  TRGSW res = trgsw_alloc_new_sample(l, Bg_bit, k, N);
  trgsw_monomial_sample(res, m, e, key);
  return res;
}

TRGSW trgsw_new_sample(Torus m, TRGSW_Key key){
  return trgsw_new_monomial_sample(m, 0, key);
}

uint64_t _debug_trgsw_decrypt_exp_sample(TRGSW c, TRGSW_Key key){
  const uint64_t N = key->trlwe_key->s[0]->N, l = key->l;
  const Torus delta = (1UL << (sizeof(Torus)*8 - 1 - key->Bg_bit));
  TorusPolynomial poly = polynomial_new_torus_polynomial(N);
  trlwe_phase(poly, c->samples[l], key->trlwe_key);
  int resIdx = -1;
  for (int j = 0; j < N; j++){
    if((poly->coeffs[j] < -delta) && (poly->coeffs[j] > delta)){
      if(resIdx != -1){
        // printf("[TRGSW Exp Decryption error] Current: %lf*x^%d - Previous: %lf*x^%d\n", torus2double(poly->coeffs[j]), j, torus2double(poly->coeffs[resIdx]), resIdx);
        // return -1;
        resIdx = -1;
        break;
      }
      resIdx = j;
    }
  }
  // if(resIdx == -1){
  //   printf("\nTRGSW error: no value\n");
  //   for (size_t i = 0; i < N; i++)
  //   {
  //     printf("%lu: %lf, ", i, torus2double(poly->coeffs[i]));
  //   }
  //   printf("\n");
  // }
  free_polynomial(poly);
  return resIdx;
}

// uint64_t _debug_trgsw_decrypt_exp_DFT_sample(TRGSW_DFT c, TRGSW_Key key){
//   const uint64_t N = key->trlwe_key->s[0]->N, k = key->trlwe_key->k, l = key->l;
//   TRLWE tmp = trlwe_new_noiseless_trivial_sample(NULL, k, N);
//   trlwe_from_DFT(tmp, c->samples[l]);
//   TorusPolynomial poly = polynomial_new_torus_polynomial(N);
//   trlwe_phase(poly, tmp, key->trlwe_key);
//   int resIdx = -1;
//   for (int j = 0; j < N; j++){
//     if((poly->coeffs[j] < - (1UL << (sizeof(Torus)*8 - 1 - key->Bg_bit))) && (poly->coeffs[j] > (1UL << (sizeof(Torus)*8 - 1 - key->Bg_bit)))){
//       if(resIdx != -1){
//         // printf("[TRGSW Exp Decryption error] Current: %lf*x^%d - Previous: %lf*x^%d\n", torus2double(poly->coeffs[j]), j, torus2double(poly->coeffs[resIdx]), resIdx);
//         return -1;
//       }
//       resIdx = j;
//     }
//   }
//   free_trlwe(tmp);
//   free_polynomial(poly);
//   return resIdx;
// }


uint64_t _debug_trgsw_decrypt_exp_DFT_sample(TRGSW_DFT c, TRGSW_Key key){
  const uint64_t N = key->trlwe_key->s[0]->N, k = key->trlwe_key->k;
  TRLWE tmp = trlwe_new_noiseless_trivial_sample(NULL, k, N);
  tmp->b->coeffs[0] = (1UL << (64 - key->Bg_bit));
  TRLWE_DFT res = trlwe_alloc_new_DFT_sample(k, N);
  trgsw_mul_trlwe_DFT(res, tmp, c);
  trlwe_from_DFT(tmp, res);
  TorusPolynomial poly = polynomial_new_torus_polynomial(N);
  trlwe_phase(poly, tmp, key->trlwe_key);
  int resIdx = -1;
  for (int j = 0; j < N; j++){
    if((poly->coeffs[j] < - (1UL << (63 - key->Bg_bit))) && (poly->coeffs[j] > (1UL << (63 - key->Bg_bit)))){
      if(resIdx != -1){
        printf("[TRGSW Exp Decryption error] Current: %lf*x^%d - Previous: %lf*x^%d\n", torus2double(poly->coeffs[j]), j, torus2double(poly->coeffs[resIdx]), resIdx);
        return -1;
      }
      resIdx = j;
    }
  }
  // if(resIdx == -1){
  //   printf("\nTRGSW error: no value\n");
  //   for (size_t i = 0; i < N; i++)
  //   {
  //     printf("%lu: %lf, ", i, torus2double(poly->coeffs[i]));
  //   }
  //   printf("\n");
  // }
  return resIdx;
}


TRGSW trgsw_new_exp_sample(int e, TRGSW_Key key){
  return trgsw_new_monomial_sample(1, e, key);
}

void trgsw_sub(TRGSW out, TRGSW in1, TRGSW in2){
  const int l = in1->l, k = in1->samples[0]->k;
  for (size_t i = 0; i < (k+1)*l; i++){
    trlwe_sub(out->samples[i], in1->samples[i], in2->samples[i]);
  }
}

void trgsw_DFT_sub(TRGSW_DFT out, TRGSW_DFT in1, TRGSW_DFT in2){
  const int l = in1->l, k = in1->samples[0]->k;
  for (size_t i = 0; i < (k+1)*l; i++){
    trlwe_DFT_sub(out->samples[i], in1->samples[i], in2->samples[i]);
  }
}

void trgsw_DFT_add(TRGSW_DFT out, TRGSW_DFT in1, TRGSW_DFT in2){
  const int l = in1->l, k = in1->samples[0]->k;
  for (size_t i = 0; i < (k+1)*l; i++){
    trlwe_DFT_add(out->samples[i], in1->samples[i], in2->samples[i]);
  }
}

void trgsw_copy(TRGSW out, TRGSW in){
  const int l = in->l, k = in->samples[0]->k;
  for (size_t i = 0; i < (k+1)*l; i++){
    trlwe_copy(out->samples[i], in->samples[i]);
  }
}

void trgsw_DFT_copy(TRGSW_DFT out, TRGSW_DFT in){
  const int l = in->l, k = in->samples[0]->k;
  for (size_t i = 0; i < (k+1)*l; i++){
    trlwe_DFT_copy(out->samples[i], in->samples[i]);
  }
}


void trgsw_add(TRGSW out, TRGSW in1, TRGSW in2){
  const int l = in1->l, k = in1->samples[0]->k;
  for (size_t i = 0; i < (k+1)*l; i++){
    trlwe_add(out->samples[i], in1->samples[i], in2->samples[i]);
  }
}

void trgsw_addto(TRGSW out, TRGSW in){
  trgsw_add(out, out, in);
}


void trgsw_mul_by_xai(TRGSW out, TRGSW in, int a){
  const int l = in->l, k = in->samples[0]->k;
  for (size_t i = 0; i < (k+1)*l; i++){
    trlwe_mul_by_xai(out->samples[i], in->samples[i], a);
  }
}

void trgsw_mul_by_xai_addto(TRGSW out, TRGSW in, int a){
  const int l = in->l, k = in->samples[0]->k;
  for (size_t i = 0; i < (k+1)*l; i++){
    trlwe_mul_by_xai_addto(out->samples[i], in->samples[i], a);
  }
}

void trgsw_mul_by_xai_minus_1(TRGSW out, TRGSW in, int a){
  const int l = in->l, k = in->samples[0]->k;
  for (size_t i = 0; i < (k+1)*l; i++){
    trlwe_mul_by_xai_minus_1(out->samples[i], in->samples[i], a);
  }
}


void trgsw_to_DFT(TRGSW_DFT out, TRGSW in){
  const int l = in->l, k = in->samples[0]->k;
  for (size_t i = 0; i < l * (k + 1); i++){
    trlwe_to_DFT(out->samples[i], in->samples[i]);
  }
}

void trgsw_from_DFT(TRGSW out, TRGSW_DFT in){
  const int l = in->l, k = in->samples[0]->k;
  for (size_t i = 0; i < l * (k + 1); i++){
    trlwe_from_DFT(out->samples[i], in->samples[i]);
  }
}

void trgsw_mul_trlwe_DFT_1(TRLWE_DFT out, TRLWE in1, TRGSW_DFT in2){
  const int N = in1->b->N, l = in2->l, k = in1->k;
  TorusPolynomial * dec_trlwe = polynomial_new_array_of_torus_polynomials(N, (k + 1) * l);
  DFT_Polynomial * dec_trlwe_DFT = polynomial_new_array_of_polynomials_DFT(N, (k + 1) * l); 

  trlwe_decompose(dec_trlwe, in1, in2->Bg_bit, in2->l);
  for (size_t i = 0; i < (k + 1) * l; i++){
    polynomial_torus_to_DFT(dec_trlwe_DFT[i], dec_trlwe[i]);
  }
  
  for (size_t i = 0; i < k; i++){
    polynomial_mul_DFT(out->a[i], dec_trlwe_DFT[0], in2->samples[0]->a[i]);
  }
  polynomial_mul_DFT(out->b, dec_trlwe_DFT[0], in2->samples[0]->b);

  for (size_t j = 1; j < (k + 1) * l; j++){
    for (size_t i = 0; i < k; i++){
      polynomial_mul_addto_DFT(out->a[i], dec_trlwe_DFT[j], in2->samples[j]->a[i]);
    }
    polynomial_mul_addto_DFT(out->b, dec_trlwe_DFT[j], in2->samples[j]->b);
  }

  free_array_of_polynomials((void **) dec_trlwe, (k + 1) * l);
  free_array_of_polynomials((void **) dec_trlwe_DFT, (k + 1) * l);
}

void trgsw_mul_trlwe_DFT(TRLWE_DFT out, TRLWE in1, TRGSW_DFT in2){
  const int N = in1->b->N, l = in2->l, k = in1->k;
  // if(in1->k > 1) return trgsw_mul_trlwe_DFT_1(out, in1, in2);
  // assert(in1->k == 1);
  assert(in1->k == in2->samples[0]->k);
  TorusPolynomial dec_poly = polynomial_new_torus_polynomial(N);
  DFT_Polynomial dec_poly_DFT = polynomial_new_DFT_polynomial(N); 

  // decomp a[0], get a[0]_0
  polynomial_decompose_i(dec_poly, in1->a[0], in2->Bg_bit, in2->l, 0);
  polynomial_torus_to_DFT(dec_poly_DFT, dec_poly);
  trlwe_DFT_mul_by_polynomial(out, in2->samples[0], dec_poly_DFT);

  // decomp a[0], get a[0]_i for i in 1 to ell - 1 
  for (size_t j = 1; j < l; j++){
    polynomial_decompose_i(dec_poly, in1->a[0], in2->Bg_bit, in2->l, j);
    polynomial_torus_to_DFT(dec_poly_DFT, dec_poly);
    trlwe_DFT_mul_addto_by_polynomial(out, in2->samples[j], dec_poly_DFT);
  }

  // decomp a[i], get a[i]_j for i in [1,k) and j in [0, ell)
  for (size_t i = 1; i < in1->k; i++){
    for (size_t j = 0; j < l; j++){
      polynomial_decompose_i(dec_poly, in1->a[i], in2->Bg_bit, in2->l, j);
      polynomial_torus_to_DFT(dec_poly_DFT, dec_poly);
      trlwe_DFT_mul_addto_by_polynomial(out, in2->samples[i*l + j], dec_poly_DFT);
    }
  }
  
  // decomp b, get b_j for j in [0, ell)
  for (size_t j = 0; j < l; j++){
    polynomial_decompose_i(dec_poly, in1->b, in2->Bg_bit, in2->l, j);
    polynomial_torus_to_DFT(dec_poly_DFT, dec_poly);
    trlwe_DFT_mul_addto_by_polynomial(out, in2->samples[k*l + j], dec_poly_DFT);
  }

  free_polynomial(dec_poly);
  free_polynomial(dec_poly_DFT);
}

void trgsw_mul_DFT(TRGSW_DFT out, TRGSW in1, TRGSW_DFT in2){
  assert(out != in2);
  const int l = in2->l, k = in1->samples[0]->k;
  for (size_t i = 0; i < (k+1)*l; i++){
    trgsw_mul_trlwe_DFT(out->samples[i], in1->samples[i], in2);
  }
}

void trgsw_mul_DFT2(TRGSW_DFT out, TRGSW_DFT in1, TRGSW_DFT in2){
  assert(out != in2);
  const int N = in1->samples[0]->b->N, l = in1->l, k = in1->samples[0]->k;
  TRLWE tmp = trlwe_alloc_new_sample(k, N);
  for (size_t i = 0; i < (k+1)*l; i++){
    trlwe_from_DFT(tmp, in1->samples[i]);
    trgsw_mul_trlwe_DFT(out->samples[i], tmp, in2);
  }
  free_trlwe(tmp);
}

/* out += in1 * in2*/
void trgsw_DFT_mul_addto_by_polynomial(TRGSW_DFT out, TRGSW_DFT in1, DFT_Polynomial in2){
  assert(out->l == in1->l);
  for (size_t i = 0; i < (in1->samples[0]->k + 1) * in1->l; i++){
    trlwe_DFT_mul_addto_by_polynomial(out->samples[i], in1->samples[i], in2);
  }
}

void trgsw_naive_mul_trlwe(TRLWE out, TRLWE in1, TRGSW in2){
  const int N = in1->b->N, l = in2->l, k = in1->k;
  TorusPolynomial * dec_trlwe = polynomial_new_array_of_torus_polynomials(N, (k + 1) * l);
  trlwe_decompose(dec_trlwe, in1, in2->Bg_bit, in2->l);
  
  for (size_t i = 0; i < k; i++){
    polynomial_naive_mul_torus(out->a[i], dec_trlwe[0], in2->samples[0]->a[i]);
  }
  polynomial_naive_mul_torus(out->b, dec_trlwe[0], in2->samples[0]->b);

  for (size_t j = 1; j < (k + 1) * l; j++){
    for (size_t i = 0; i < k; i++){
      polynomial_naive_mul_addto_torus(out->a[i], dec_trlwe[j], in2->samples[j]->a[i]);
    }
    polynomial_naive_mul_addto_torus(out->b, dec_trlwe[j], in2->samples[j]->b);
  }

  free_array_of_polynomials((void **) dec_trlwe, (k + 1) * l);
}

void trgsw_naive_mul(TRGSW out, TRGSW in1, TRGSW in2){
  const int l = in2->l, k = in1->samples[0]->k;
  for (size_t i = 0; i < (k+1)*l; i++){
    trgsw_naive_mul_trlwe(out->samples[i], in1->samples[i], in2);
  }
}

void trgsw_ks_b_to_a(TRGSW c, TRLWE_KS_Key * ksk){
  for (size_t i = 0; i < c->l; i++){
    trlwe_priv_keyswitch_2(c->samples[i], c->samples[c->l + i], ksk);
  }
}