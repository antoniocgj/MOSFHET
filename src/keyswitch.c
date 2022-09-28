#include "mosfhet.h"

TRLWE_KS_Key trlwe_new_RL_key(TRLWE_Key key, int t, int base_bit){
  assert(key->k == 1);
  TRLWE_Key key2 = trlwe_alloc_key(key->s[0]->N, key->k, key->sigma);
  polynomial_naive_mul_torus(key2->s[0], key->s[0], key->s[0]);
  TRLWE_KS_Key res = trlwe_new_KS_key(key, key2, t, base_bit);
  free_trlwe_key(key2);
  return res;
}

TRLWE_KS_Key trlwe_new_KS_key(TRLWE_Key out_key, TRLWE_Key in_key, int t, int base_bit){
  const int bit_size = sizeof(Torus)*8, N_out = out_key->s[0]->N, N_in = in_key->s[0]->N;
  TRLWE_KS_Key res;
  res = (TRLWE_KS_Key) safe_malloc(sizeof(*res));
  res->base_bit = base_bit;
  res->k = in_key->k;
  res->t = t;
  TorusPolynomial dec_poly = polynomial_new_torus_polynomial(N_in);
  TRLWE tmp = trlwe_alloc_new_sample(out_key->k, N_out);

  res->s = (TRLWE_DFT **) safe_malloc(sizeof(TRLWE_DFT*) * in_key->k);
  for (size_t i = 0; i < in_key->k; i++){
    res->s[i] = (TRLWE_DFT *) safe_malloc(sizeof(TRLWE_DFT) * t);
    for (size_t j = 0; j < t; j++){
      for (size_t i2 = 0; i2 < N_in; i2++){
        dec_poly->coeffs[i2] = in_key->s[i]->coeffs[i2] * (1UL << (bit_size - (j + 1) * base_bit));
      }
      trlwe_sample(tmp, dec_poly, out_key);
      res->s[i][j] = trlwe_alloc_new_DFT_sample(out_key->k, N_out);
      trlwe_to_DFT(res->s[i][j], tmp);
    }
  }
  free_polynomial(dec_poly);
  free_trlwe(tmp);
  return res;
}

TRLWE_KS_Key * trlwe_new_priv_KS_key(TRLWE_Key out_key, TRLWE_Key in_key, int t, int base_bit){
  assert(out_key->k == 1);
  TRLWE_Key tmp_key = trlwe_alloc_key(out_key->s[0]->N, out_key->k, out_key->sigma);
  polynomial_negate_torus_polynomial(tmp_key->s[0], out_key->s[0]);
  polynomial_mul_torus(tmp_key->s[0], tmp_key->s[0], in_key->s[0]);
  TRLWE_KS_Key * res = safe_malloc(sizeof(TRLWE_KS_Key)*2);
  res[0] = trlwe_new_KS_key(out_key, tmp_key, t, base_bit);
  polynomial_negate_torus_polynomial(tmp_key->s[0], out_key->s[0]);
  res[1] = trlwe_new_KS_key(out_key, tmp_key, t, base_bit);
  free_trlwe_key(tmp_key);
  return res;
}

void trlwe_priv_keyswitch_2(TRLWE out, TRLWE in, TRLWE_KS_Key * ks_key){
  TRLWE tmp = trlwe_alloc_new_sample(in->k, in->b->N);
  polynomial_negate_torus_polynomial(tmp->a[0], in->b);
  memset(tmp->b->coeffs, 0, sizeof(Torus)*tmp->b->N);
  trlwe_keyswitch(tmp, tmp, ks_key[1]);
  polynomial_copy_torus_polynomial(out->a[0], in->a[0]);
  memset(out->b->coeffs, 0, sizeof(Torus)*out->b->N);
  trlwe_keyswitch(out, out, ks_key[0]);
  trlwe_addto(out, tmp);
  free_trlwe(tmp);
}

TRLWE_KS_Key trlwe_new_full_packing_KS_key(TRLWE_Key out_key, TLWE_Key in_key, int t, int base_bit){
  TRLWE_Key tmp = trlwe_alloc_key(1, in_key->n, in_key->sigma); // N=1, k = n 
  for (size_t i = 0; i < in_key->n; i++){
    tmp->s[i]->coeffs[0] = in_key->s[i];
  }
  TRLWE_KS_Key res = trlwe_new_KS_key(out_key, tmp, t, base_bit);
  free_trlwe_key(tmp);
  return res;
}

void free_trlwe_ks_key(TRLWE_KS_Key key){
  const int t = key->t, k = key->k;
  for (size_t i = 0; i < k; i++){
    for (size_t j = 0; j < t; j++){
      free_trlwe(key->s[i][j]);
    }
    free(key->s[i]);
  }
  free(key->s);
  free(key);
}


void trlwe_save_KS_key(FILE * fd, TRLWE_KS_Key key){
  const int t = key->t, k_in = key->k, N = key->s[0][0]->b->N, k = key->s[0][0]->k;
  fwrite(&key->base_bit, sizeof(int), 1, fd);
  fwrite(&t, sizeof(int), 1, fd);
  fwrite(&k_in, sizeof(int), 1, fd);
  fwrite(&k, sizeof(int), 1, fd);
  fwrite(&N, sizeof(int), 1, fd);
  for (size_t i = 0; i < k; i++){
    for (size_t j = 0; j < t; j++){
      trlwe_save_DFT_sample(fd, key->s[i][j]);
    }
  }
}

TRLWE_KS_Key trlwe_load_new_KS_key(FILE * fd){
  int base_bit, t, k_in, N, k;
  fread(&base_bit, sizeof(int), 1, fd);
  fread(&t, sizeof(int), 1, fd);
  fread(&k_in, sizeof(int), 1, fd);
  fread(&k, sizeof(int), 1, fd);
  fread(&N, sizeof(int), 1, fd);

  TRLWE_KS_Key res;
  res = (TRLWE_KS_Key) safe_malloc(sizeof(*res));
  res->base_bit = base_bit;
  res->t = t;
  res->k = k_in;


  res->s = (TRLWE_DFT **) safe_malloc(sizeof(TRLWE_DFT*) * k_in);
  for (size_t i = 0; i < k_in; i++){
    res->s[i] = (TRLWE_DFT *) safe_malloc(sizeof(TRLWE_DFT) * t);
    for (size_t j = 0; j < t; j++){
      res->s[i][j] = trlwe_load_new_DFT_sample(fd, k, N);
    }
  }
  return res;
}


void trlwe_keyswitch(TRLWE out, TRLWE in, TRLWE_KS_Key ks_key){
  const int N = out->b->N;
  assert(out->k == ks_key->s[0][0]->k);
  assert(out->b->N == ks_key->s[0][0]->b->N);
  

  TorusPolynomial dec_in_a = polynomial_new_torus_polynomial(N);
  DFT_Polynomial tmp = polynomial_new_DFT_polynomial(N);
  TRLWE_DFT acc = trlwe_alloc_new_DFT_sample(out->k, N);
  TRLWE as = trlwe_alloc_new_sample(out->k, N);

  for (size_t i = 0; i < out->k; i++){
    for (size_t j = 0; j < N; j++) acc->a[i]->coeffs[j] = 0.;
  }
  for (size_t j = 0; j < N; j++) acc->b->coeffs[j] = 0.;
  
  for (size_t i = 0; i < in->k; i++){
    for (size_t j = 0; j < ks_key->t; j++){
      polynomial_decompose_i(dec_in_a, in->a[i], ks_key->base_bit, ks_key->t, j);
      polynomial_torus_to_DFT(tmp, dec_in_a);
      trlwe_DFT_mul_addto_by_polynomial(acc, ks_key->s[i][j], tmp);   
    }
  }
  trlwe_from_DFT(as, acc);
  trlwe_noiseless_trivial_sample(out, in->b);
  trlwe_subto(out, as);
  
  free_trlwe(acc);
  free_trlwe(as);
  free_polynomial(tmp);
  free_polynomial(dec_in_a);
}

void trlwe_full_packing_keyswitch(TRLWE out, TLWE * in, uint64_t size, TRLWE_KS_Key ks_key){
  const int N = out->b->N;
  assert(out->k == ks_key->s[0][0]->k);
  assert(out->b->N == ks_key->s[0][0]->b->N);

  TorusPolynomial dec_in_a = polynomial_new_torus_polynomial(N);
  TorusPolynomial a_i = polynomial_new_torus_polynomial(N);
  DFT_Polynomial tmp = polynomial_new_DFT_polynomial(N);
  TRLWE_DFT acc = trlwe_alloc_new_DFT_sample(out->k, N);

  memset(a_i->coeffs, 0, sizeof(Torus)*N);
  for (size_t i = 0; i < out->k; i++){
    memset(acc->a[i]->coeffs, 0, sizeof(Torus)*N);
  }
  memset(acc->b->coeffs, 0, sizeof(Torus)*N);
  
  for (size_t i = 0; i < in[0]->n; i++){
    for (size_t j = 0; j < size; j++){
      a_i->coeffs[j] = in[j]->a[i];
    }
    for (size_t j = 0; j < ks_key->t; j++){
      polynomial_decompose_i(dec_in_a, a_i, ks_key->base_bit, ks_key->t, j);
      polynomial_torus_to_DFT(tmp, dec_in_a);
      trlwe_DFT_mul_addto_by_polynomial(acc, ks_key->s[i][j], tmp);   
    }
  }
  trlwe_from_DFT(out, acc);
  trlwe_negate(out, out);
  for (size_t j = 0; j < size; j++){
    out->b->coeffs[j] += in[j]->b;
  }
  
  free_trlwe(acc);
  free_polynomial(tmp);
  free_polynomial(dec_in_a);
}


/* Packing Key Switching */
#ifdef USE_COMPRESSED_TRLWE
#define _MACRO_trlwe_new_sample trlwe_new_compressed_sample
#define _MACRO_trlwe_save_sample trlwe_save_compressed_sample
#define _MACRO_trlwe_load_new_sample trlwe_load_new_compressed_sample
#define _MACRO_trlwe_subto trlwe_compressed_subto
#else
#define _MACRO_trlwe_new_sample trlwe_new_sample
#define _MACRO_trlwe_save_sample trlwe_save_sample
#define _MACRO_trlwe_load_new_sample trlwe_load_new_sample
#define _MACRO_trlwe_subto trlwe_subto
#endif


LUT_Packing_KS_Key trlwe_new_packing_KS_key(TRLWE_Key out_key, TLWE_Key in_key, int t, int base_bit, int torus_base){
  const int base = 1 << base_bit, bit_size = sizeof(Torus)*8, N = out_key->s[0]->N;
  LUT_Packing_KS_Key res;
  res = (LUT_Packing_KS_Key) safe_malloc(sizeof(*res));
  res->base_bit = base_bit;
  res->t = t;
  res->n = in_key->n;
  res->torus_base = torus_base;


  res->s = (TRLWE ****) safe_malloc(sizeof(TRLWE****) * in_key->n);
  for (size_t i = 0; i < in_key->n; i++){
    res->s[i] = (TRLWE ***) safe_malloc(sizeof(TRLWE***) * torus_base);
    for (size_t e = 0; e < torus_base; e++){
      res->s[i][e] = (TRLWE **) safe_malloc(sizeof(TRLWE**) * t);
      for (size_t j = 0; j < t; j++){
        res->s[i][e][j] = (TRLWE*) safe_malloc(sizeof(TRLWE*) * (base - 1));
        for (size_t k = 0; k < base - 1; k++){
          const Torus dec_key = in_key->s[i] * (k + 1) * (1UL << (bit_size - (j + 1) * base_bit));
          res->s[i][e][j][k] = _MACRO_trlwe_new_sample(0, out_key);
          for (size_t q = e*(N/torus_base); q < (e+1)*(N/torus_base); q++) res->s[i][e][j][k]->b->coeffs[q] += dec_key;
        }
      }
    }
  }
  return res;
}

void trlwe_save_packing_KS_key(FILE * fd, LUT_Packing_KS_Key key){
  const int base = 1 << key->base_bit, t = key->t, torus_base = key->torus_base, n = key->n, N = key->s[0][0][0][0]->b->N, k = key->s[0][0][0][0]->k;
  fwrite(&key->base_bit, sizeof(int), 1, fd);
  fwrite(&t, sizeof(int), 1, fd);
  fwrite(&torus_base, sizeof(int), 1, fd);
  fwrite(&n, sizeof(int), 1, fd);
  fwrite(&k, sizeof(int), 1, fd);
  fwrite(&N, sizeof(int), 1, fd);
  for (size_t i = 0; i < n; i++){
    for (size_t e = 0; e < torus_base; e++){
      for (size_t j = 0; j < t; j++){
        for (size_t l = 0; l < base - 1; l++){
          _MACRO_trlwe_save_sample(fd, key->s[i][e][j][l]);
        }
      }
    }
  }
}

LUT_Packing_KS_Key trlwe_load_new_packing_KS_key(FILE * fd){
  int base_bit, t, torus_base, n, N, k;
  fread(&base_bit, sizeof(int), 1, fd);
  fread(&t, sizeof(int), 1, fd);
  fread(&torus_base, sizeof(int), 1, fd);
  fread(&n, sizeof(int), 1, fd);
  fread(&k, sizeof(int), 1, fd);
  fread(&N, sizeof(int), 1, fd);
  const int base = 1 << base_bit;

  LUT_Packing_KS_Key res;
  res = (LUT_Packing_KS_Key) safe_malloc(sizeof(*res));
  res->base_bit = base_bit;
  res->t = t;
  res->n = n;
  res->torus_base = torus_base;


  res->s = (TRLWE ****) safe_malloc(sizeof(TRLWE****) * n);
  for (size_t i = 0; i < n; i++){
    res->s[i] = (TRLWE ***) safe_malloc(sizeof(TRLWE***) * torus_base);
    for (size_t e = 0; e < torus_base; e++){
      res->s[i][e] = (TRLWE **) safe_malloc(sizeof(TRLWE**) * t);
      for (size_t j = 0; j < t; j++){
        res->s[i][e][j] = (TRLWE*) safe_malloc(sizeof(TRLWE*) * (base - 1));
        for (size_t l = 0; l < base - 1; l++){
          res->s[i][e][j][l] = _MACRO_trlwe_load_new_sample(fd, k, N);
        }
      }
    }
  }
  return res;
}

void free_trlwe_packing_ks_key(LUT_Packing_KS_Key key){
  const int base = 1 << key->base_bit, t = key->t, torus_base = key->torus_base, n = key->n;
  for (size_t i = 0; i < n; i++){
    for (size_t e = 0; e < torus_base; e++){
      for (size_t j = 0; j < t; j++){
        for (size_t k = 0; k <  base - 1; k++){
          free_trlwe(key->s[i][e][j][k]);
        }
        free(key->s[i][e][j]);
      }
      free(key->s[i][e]);
    }
    free(key->s[i]);
  }
  free(key->s);
  free(key);
}

void trlwe_packing_keyswitch(TRLWE out, TLWE * in, LUT_Packing_KS_Key ks_key){
  const int bit_size = sizeof(Torus)*8, N = out->b->N, torus_base = ks_key->torus_base;
  const Torus prec_offset = 1UL << (bit_size - (1 + ks_key->base_bit * ks_key->t));
  const Torus mask = (1UL << ks_key->base_bit) - 1;
  assert(out->k == ks_key->s[0][0][0][0]->k);
  assert(out->b->N == ks_key->s[0][0][0][0]->b->N);

  for (size_t i = 0; i < N; i++) {
    out->a[0]->coeffs[i] = 0;
    out->b->coeffs[i] = in[i/(N/torus_base)]->b;
  }

  for (size_t i = 0; i < in[0]->n;i++){
    for (size_t e = 0; e < torus_base; e++){
      const Torus aibar = in[e]->a[i]+prec_offset;
      for (size_t j = 0; j < ks_key->t; j++){
        const Torus aij = (aibar>>(bit_size-(j+1)*ks_key->base_bit)) & mask;
        if(aij != 0) _MACRO_trlwe_subto(out, ks_key->s[i][e][j][aij - 1]);
      }
    }     
  }
}

/* TLWE(m) -> TRLWE(m*X^0) Key Switching */

Generic_KS_Key trlwe_new_packing1_KS_key(TRLWE_Key out_key, TLWE_Key in_key, int t, int base_bit){
  const int base = 1 << base_bit, bit_size = sizeof(Torus)*8;
  Generic_KS_Key res;
  res = (Generic_KS_Key) safe_malloc(sizeof(*res));
  res->base_bit = base_bit;
  res->t = t;
  res->n = in_key->n;
  res->include_b = 0;

  res->s = (TRLWE ***) safe_malloc(sizeof(TRLWE****) * in_key->n);
  for (size_t i = 0; i < in_key->n; i++){
    res->s[i] = (TRLWE **) safe_malloc(sizeof(TRLWE**) * t);
    for (size_t j = 0; j < t; j++){
      res->s[i][j] = (TRLWE*) safe_malloc(sizeof(TRLWE*) * (base - 1));
      for (size_t k = 0; k < base - 1; k++){
        const Torus dec_key = in_key->s[i] * (k + 1) * (1UL << (bit_size - (j + 1) * base_bit));
        res->s[i][j][k] = _MACRO_trlwe_new_sample(0, out_key);
        res->s[i][j][k]->b->coeffs[0] += dec_key;
      }
    }
  }
  return res;
}


void free_trlwe_generic_ks_key(Generic_KS_Key key){
  const int base = 1 << key->base_bit, t = key->t, n = key->n;
  for (size_t i = 0; i < n + key->include_b; i++){
    for (size_t j = 0; j < t; j++){
      for (size_t k = 0; k <  base - 1; k++){
        free_trlwe(key->s[i][j][k]);
      }
      free(key->s[i][j]);
    }
    free(key->s[i]);
  }
  free(key->s);
  free(key);
}


void trlwe_save_generic_ks_key(FILE * fd, Generic_KS_Key key){
  const int base = 1 << key->base_bit, t = key->t, n = key->n, N = key->s[0][0][0]->b->N, k = key->s[0][0][0]->k, include_b = key->include_b;
  fwrite(&key->base_bit, sizeof(int), 1, fd);
  fwrite(&t, sizeof(int), 1, fd);
  fwrite(&n, sizeof(int), 1, fd);
  fwrite(&k, sizeof(int), 1, fd);
  fwrite(&N, sizeof(int), 1, fd);
  fwrite(&include_b, sizeof(int), 1, fd);
  for (size_t i = 0; i < n + include_b; i++){
    for (size_t j = 0; j < t; j++){
      for (size_t l = 0; l < base - 1; l++){
        _MACRO_trlwe_save_sample(fd, key->s[i][j][l]);
      }
    }
  }
}

Generic_KS_Key trlwe_load_new_generic_ks_key(FILE * fd){
  int base_bit, t, n, N, k, include_b;
  fread(&base_bit, sizeof(int), 1, fd);
  fread(&t, sizeof(int), 1, fd);
  fread(&n, sizeof(int), 1, fd);
  fread(&k, sizeof(int), 1, fd);
  fread(&N, sizeof(int), 1, fd);
  fread(&include_b, sizeof(int), 1, fd);

  const int base = 1 << base_bit;

  Generic_KS_Key res;
  res = (Generic_KS_Key) safe_malloc(sizeof(*res));
  res->base_bit = base_bit;
  res->t = t;
  res->n = n;
  res->include_b = include_b;

  res->s = (TRLWE ***) safe_malloc(sizeof(TRLWE**) * n);
  for (size_t i = 0; i < n + include_b; i++){
    res->s[i] = (TRLWE **) safe_malloc(sizeof(TRLWE*) * t);
      for (size_t j = 0; j < t; j++){
        res->s[i][j] = (TRLWE*) safe_malloc(sizeof(TRLWE) * (base - 1));
        for (size_t l = 0; l < base - 1; l++){
          res->s[i][j][l] = _MACRO_trlwe_load_new_sample(fd, k, N);
        }
      }
  }
  return res;
}


void trlwe_packing1_keyswitch(TRLWE out, TLWE in, Generic_KS_Key ks_key){
  const int bit_size = sizeof(Torus)*8;
  const Torus prec_offset = 1UL << (bit_size - (1 + ks_key->base_bit * ks_key->t));
  const Torus mask = (1UL << ks_key->base_bit) - 1;
  assert(out->k == ks_key->s[0][0][0]->k);
  assert(out->b->N == ks_key->s[0][0][0]->b->N);

  trlwe_noiseless_trivial_sample(out, NULL);
  out->b->coeffs[0] = in->b;

  for (size_t i = 0; i < in->n; i++){
    const Torus aibar = in->a[i]+prec_offset;
    for (size_t j = 0; j < ks_key->t; j++){
      const Torus aij = (aibar>>(bit_size-(j+1)*ks_key->base_bit)) & mask;
      if(aij != 0) _MACRO_trlwe_subto(out, ks_key->s[i][j][aij - 1]);
    }     
  }
}

TRLWE_KS_Key * trlwe_new_packing1_KS_key_CDKS21(TRLWE_Key out_key, TLWE_Key in_key, int t, int base_bit){
  const int N = out_key->s[0]->N, log_N = log2(N);
  TRLWE_Key key2 = trlwe_alloc_key(N, 1, in_key->sigma);
  TorusPolynomial tmp = polynomial_new_torus_polynomial(N);
  TorusPolynomial tmp2 = polynomial_new_torus_polynomial(N);
  for (size_t i = 0; i < in_key->n; i++) tmp->coeffs[i] = in_key->s[i];
  for (size_t i = in_key->n; i < N; i++) tmp->coeffs[i] = 0;
  TRLWE_KS_Key * res = (TRLWE_KS_Key * ) safe_malloc(sizeof(TRLWE_KS_Key) * log_N);

  for (size_t j = 0; j < log_N; j++){
    polynomial_permute(tmp2, tmp, (1 << (log_N - j)) + 1);
    for (size_t i = 0; i < N; i++){
      key2->s[0]->coeffs[i] = tmp2->coeffs[i];
    }
    res[j] = trlwe_new_KS_key(out_key, key2, t, base_bit);
  }
  
  free_trlwe_key(key2);
  free_polynomial(tmp);
  free_polynomial(tmp2);
  return res;
}

TRLWE_KS_Key * trlwe_new_automorphism_KS_keyset(TRLWE_Key key, bool skip_even, int t, int base_bit){
  const int N = key->s[0]->N;
  TRLWE_Key key2 = trlwe_alloc_key(N, 1, key->sigma);
  TRLWE_KS_Key * res = (TRLWE_KS_Key * ) safe_malloc(sizeof(TRLWE_KS_Key) * N);
  for (size_t i = 0, j = 0; i < (N<<1); i++){
    if(skip_even && !(i&1)) continue;
    polynomial_permute(key2->s[0], key->s[0], i);
    res[j++] = trlwe_new_KS_key(key, key2, t, base_bit);
  }  
  free_trlwe_key(key2);
  return res;
}

void trlwe_packing1_keyswitch_CDKS21(TRLWE out, TLWE in, TRLWE_KS_Key * ks_key){
  const int k = out->k, N = out->b->N;
  assert(out->k == ks_key[0]->k);
  assert(out->b->N == ks_key[0]->s[0][0]->b->N);
  TRLWE tmp = trlwe_alloc_new_sample(k, N);
  trlwe_noiseless_trivial_sample(out, 0);
  // T^n -> T_N[X]
  for (size_t i = 1; i < N; i++) out->a[0]->coeffs[N - i] = -in->a[i];
  out->a[0]->coeffs[0] = in->a[0];
  // memset(out->b->coeffs, 0, sizeof(Torus)*N);
  out->b->coeffs[0] = in->b;
  // Trace
  for (size_t i = 1, j = 0; i < N; i<<=1, j++){
    const uint64_t gen = (N>>j) + 1;
    polynomial_permute(tmp->a[0], out->a[0], gen);
    polynomial_permute(tmp->b, out->b, gen);
    trlwe_keyswitch(tmp, tmp, ks_key[j]);
    trlwe_addto(out, tmp);
  }  
  free_trlwe(tmp);
}

/* Private Key Switching TLWE(M) -> TRLWE(m*-s) */
Generic_KS_Key trlwe_new_priv_SK_KS_key(TRLWE_Key out_key, TLWE_Key in_key, int t, int base_bit){
  const int base = 1 << base_bit, bit_size = sizeof(Torus)*8;
  assert(out_key->k == 1 /* new TRLWE definition */);
  Generic_KS_Key res;
  res = (Generic_KS_Key) safe_malloc(sizeof(*res));
  res->base_bit = base_bit;
  res->t = t;
  res->n = in_key->n;
  res->include_b = 1;

  res->s = (TRLWE ***) safe_malloc(sizeof(TRLWE****) * (in_key->n + 1));
  for (size_t i = 0; i < in_key->n + 1; i++){
    const Torus s_i = i < in_key->n ? in_key->s[i] : -1;
    res->s[i] = (TRLWE **) safe_malloc(sizeof(TRLWE**) * t);
    for (size_t j = 0; j < t; j++){
      res->s[i][j] = (TRLWE*) safe_malloc(sizeof(TRLWE*) * (base - 1));
      for (size_t k = 0; k < base - 1; k++){
        const Torus dec_key = s_i * (k + 1) * (1UL << (bit_size - (j + 1) * base_bit));
        res->s[i][j][k] = _MACRO_trlwe_new_sample(0, out_key);
        for (size_t e = 0; e < out_key->s[0]->N; e++){
          res->s[i][j][k]->b->coeffs[e] += (-out_key->s[0]->coeffs[e])*dec_key;
        }
      }
    }
  }
  return res;
}

void trlwe_priv_keyswitch(TRLWE out, TLWE in, Generic_KS_Key ks_key){
  const int bit_size = sizeof(Torus)*8;
  const Torus prec_offset = 1UL << (bit_size - (1 + ks_key->base_bit * ks_key->t));
  const Torus mask = (1UL << ks_key->base_bit) - 1;
  assert(out->k == ks_key->s[0][0][0]->k);
  assert(out->b->N == ks_key->s[0][0][0]->b->N);

  trlwe_noiseless_trivial_sample(out, NULL);

  for (size_t i = 0; i < in->n + 1; i++){
    const Torus a_i = i < in->n ? in->a[i] : in->b;
    const Torus aibar = a_i+prec_offset;
    for (size_t j = 0; j < ks_key->t; j++){
      const Torus aij = (aibar>>(bit_size-(j+1)*ks_key->base_bit)) & mask;
      if(aij != 0) _MACRO_trlwe_subto(out, ks_key->s[i][j][aij - 1]);
    }     
  }
}
