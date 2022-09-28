#include "mosfhet.h"

TLWE tlwe_alloc_sample(int n){
  TLWE res;
  res = (TLWE) safe_malloc(sizeof(*res));
  res->a = (Torus *) safe_aligned_malloc(sizeof(Torus) * n);
  res->n = n;
  return res;
}

TLWE * tlwe_alloc_sample_array(int count, int n){
  TLWE * res = (TLWE *) safe_malloc(sizeof(TLWE) * count);
  for (size_t i = 0; i < count; i++){
    res[i] = tlwe_alloc_sample(n);
  }
  return res;
}

TLWE tlwe_new_noiseless_trivial_sample(Torus m, int n){
  TLWE res = tlwe_alloc_sample(n);
  memset(res->a, 0, sizeof(Torus) * n);
  res->b = m;
  return res;
}

void tlwe_noiseless_trivial_sample(TLWE out, Torus m){
  memset(out->a, 0, sizeof(Torus) * out->n);
  out->b = m;
}

void free_tlwe_array(TLWE * p, int count){
  for (size_t i = 0; i < count; i++){
    free_tlwe(p[i]);
  }  
  free(p);
}

void free_tlwe(TLWE p){
  free(p->a);
  free(p);
}

void tlwe_save_sample(FILE * fd, TLWE c){
  fwrite(c->a, sizeof(Torus), c->n, fd);
  fwrite(&c->b, sizeof(Torus), 1, fd);
}

TLWE tlwe_load_new_sample(FILE * fd, int n){
  TLWE res = tlwe_alloc_sample(n);
  fread(res->a, sizeof(Torus), n, fd);
  fread(&res->b, sizeof(Torus), 1, fd);
  return res;
}

void tlwe_load_sample(FILE * fd, TLWE c){
  fread(c->a, sizeof(Torus), c->n, fd);
  fread(&c->b, sizeof(Torus), 1, fd);
}

TLWE_Key tlwe_alloc_key(int n, double sigma){
  TLWE_Key res;
  res = (TLWE_Key) safe_malloc(sizeof(*res));
  res->n = n;
  res->sigma = sigma;
  res->s = (Integer *) safe_aligned_malloc(sizeof(Integer)*n);
  return res;
}

// bound must be a power of 2
TLWE_Key tlwe_new_bounded_key(int n, uint64_t bound, double sigma){
  TLWE_Key res = tlwe_alloc_key(n, sigma);
  generate_random_bytes(n*sizeof(Integer), (uint8_t *) res->s);
  for (size_t i = 0; i < n; i++){
    res->s[i] &= (bound - 1);
    res->s[i] -= (bound >> 1) - 1;
  }
  return res;
}

TLWE_Key tlwe_new_binary_key(int n, double sigma){
  return tlwe_new_bounded_key(n, 2, sigma);
}


void tlwe_save_key(FILE * fd, TLWE_Key key){
  fwrite(&key->n, sizeof(int), 1, fd);
  fwrite(&key->sigma, sizeof(double), 1, fd);
  fwrite(key->s, sizeof(Torus), key->n, fd);
}

TLWE_Key tlwe_load_new_key(FILE * fd){
  int n;
  double sigma;
  fread(&n, sizeof(int), 1, fd);
  fread(&sigma, sizeof(double), 1, fd);
  TLWE_Key key = tlwe_alloc_key(n, sigma);
  fread(key->s, sizeof(Torus), key->n, fd);
  return key;
}

void free_tlwe_key(TLWE_Key key){
  free(key->s);
  free(key);
}

void tlwe_sample(TLWE out, Torus m, TLWE_Key key){
  const int byte_size = sizeof(Torus) * key->n;
  generate_random_bytes(byte_size, (uint8_t *) out->a);
  out->b = m;
  // internal product
  for (size_t i = 0; i < key->n; i++){
    out->b += key->s[i] * out->a[i];
  }
  out->b += double2torus(generate_normal_random(key->sigma));
}

void tlwe_copy(TLWE out, TLWE in){
  memcpy(out->a, in->a, sizeof(Torus)*in->n);
  out->b = in->b;
}

TLWE tlwe_new_sample(Torus m, TLWE_Key key){
  TLWE res = tlwe_alloc_sample(key->n);
  const int byte_size = sizeof(Torus) * key->n;
  generate_random_bytes(byte_size, (uint8_t *) res->a);
  res->b = m;
  // internal product
  for (size_t i = 0; i < key->n; i++){
    res->b += key->s[i] * res->a[i];
  }
  res->b += double2torus(generate_normal_random(key->sigma));
  return res;
}

Torus tlwe_phase(TLWE c, TLWE_Key key){
  Torus sa = 0;
  for (size_t i = 0; i < key->n; i++){
    sa += key->s[i] * c->a[i];
  }
  return c->b - sa;
}

void tlwe_add(TLWE out, TLWE in1, TLWE in2){
  for (size_t i = 0; i < in1->n; i++){
    out->a[i] = in1->a[i] + in2->a[i];
  }
  out->b = in1->b + in2->b;
}

void tlwe_addto(TLWE out, TLWE in){
  tlwe_add(out, out, in);
}

void tlwe_sub(TLWE out, TLWE in1, TLWE in2){
  for (size_t i = 0; i < in1->n; i++){
    out->a[i] = in1->a[i] - in2->a[i];
  }
  out->b = in1->b - in2->b;
}

void tlwe_negate(TLWE out, TLWE in){
  for (size_t i = 0; i < in->n; i++){
    out->a[i] = -in->a[i];
  }
  out->b = -in->b;
}

void tlwe_subto(TLWE out, TLWE in){
  tlwe_sub(out, out, in);
}

TLWE_KS_Key tlwe_new_KS_key(TLWE_Key out_key, TLWE_Key in_key, int t, int base_bit){
  const int base = 1 << base_bit, bit_size = sizeof(Torus)*8;
  TLWE_KS_Key res;
  res = (TLWE_KS_Key) safe_malloc(sizeof(*res));
  res->base_bit = base_bit;
  res->t = t;
  res->n = in_key->n;

  res->s = (TLWE ***) safe_malloc(sizeof(TLWE**) * in_key->n);
  for (size_t i = 0; i < in_key->n; i++){
    res->s[i] = (TLWE **) safe_malloc(sizeof(TLWE*) * t);
    for (size_t j = 0; j < t; j++){
      res->s[i][j] = (TLWE*) safe_malloc(sizeof(TLWE) * (base - 1));
      for (size_t k = 0; k < base - 1; k++){
        res->s[i][j][k] = tlwe_new_sample(in_key->s[i] * (k + 1) * (1ULL << (bit_size - (j + 1) * base_bit)), out_key);
      }
    }
  }
  return res;
}

void free_tlwe_ks_key(TLWE_KS_Key key){
  const int base = 1 << key->base_bit, t = key->t, n = key->n;
  for (size_t i = 0; i < n; i++){
    for (size_t j = 0; j < t; j++){
      for (size_t k = 0; k < base - 1; k++){
        free_tlwe(key->s[i][j][k]);
      }
      free(key->s[i][j]);
    }
    free(key->s[i]);
  }
  free(key->s);
  free(key);
}

TLWE_KS_Key tlwe_load_new_KS_key(FILE * fd){
  int n, t, base_bit, n_outkey;
  fread(&n, 1, sizeof(int), fd);
  fread(&t, 1, sizeof(int), fd);
  fread(&base_bit, 1, sizeof(int), fd);
  fread(&n_outkey, 1, sizeof(int), fd);
  const int base = 1 << base_bit;

  TLWE_KS_Key res;
  res = (TLWE_KS_Key) safe_malloc(sizeof(*res));
  res->base_bit = base_bit;
  res->t = t;
  res->n = n;

  res->s = (TLWE ***) safe_malloc(sizeof(TLWE**) * n);
  for (size_t i = 0; i < n; i++){
    res->s[i] = (TLWE **) safe_malloc(sizeof(TLWE*) * t);
    for (size_t j = 0; j < t; j++){
      res->s[i][j] = (TLWE*) safe_malloc(sizeof(TLWE) * (base - 1));
      for (size_t k = 0; k < base - 1; k++){
        res->s[i][j][k] = tlwe_load_new_sample(fd, n_outkey);
      }
    }
  }
  return res;
}

void tlwe_save_KS_key(FILE * fd, TLWE_KS_Key key){
  fwrite(&key->n, 1, sizeof(int), fd);
  fwrite(&key->t, 1, sizeof(int), fd);
  fwrite(&key->base_bit, 1, sizeof(int), fd);
  fwrite(&key->s[0][0][0]->n, 1, sizeof(int), fd);
  
  for (size_t i = 0; i < key->n; i++){
    for (size_t j = 0; j < key->t; j++){
      for (size_t k = 0; k < (1 << key->base_bit) - 1; k++){
        tlwe_save_sample(fd, key->s[i][j][k]);
      }
    }
  }
}

void tlwe_keyswitch(TLWE out, TLWE in, TLWE_KS_Key ks_key){
  const int bit_size = sizeof(Torus)*8;
  const Torus prec_offset = 1ULL << (bit_size - (1 + ks_key->base_bit * ks_key->t));
  const Torus mask = (1ULL << ks_key->base_bit) - 1;
  assert(out->n == ks_key->s[0][0][0]->n);

  tlwe_noiseless_trivial_sample(out, in->b);
  for (size_t i = 0; i < in->n; i++) {
    const Torus ai = in->a[i] + prec_offset;
    for (size_t j = 0; j < ks_key->t; j++) {
      const Torus aij = (ai >> (bit_size - (j + 1) * ks_key->base_bit)) & mask;
      if (aij != 0) tlwe_subto(out, ks_key->s[i][j][aij - 1]);
    }
  }
}

void tlwe_mul(TLWE out, TLWE in1, TLWE in2, int precision, Generic_KS_Key ksk, TRLWE_KS_Key rlk){
  const int N = ksk->s[0][0][0]->b->N;
  TRLWE tmp1 = trlwe_alloc_new_sample(1, N);
  TRLWE tmp2 = trlwe_alloc_new_sample(1, N);
  trlwe_packing1_keyswitch(tmp1, in1, ksk);
  trlwe_packing1_keyswitch(tmp2, in2, ksk);
  trlwe_tensor_prod_FFT(tmp1, tmp1, tmp2, precision, rlk);
  trlwe_extract_tlwe(out, tmp1, 0);
  free_trlwe(tmp1);
  free_trlwe(tmp2);
}