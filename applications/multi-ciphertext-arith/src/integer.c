#include "ufhe.h"
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

ufhe_integer ufhe_new_integer(int precision, bool _signed, ufhe_context ctx){
  const int d = precision/ctx->log_torus_base + ((precision%ctx->log_torus_base) != 0);
  ufhe_integer res;
  res = (ufhe_integer) safe_malloc(sizeof(*res));
  res->digits = tlwe_alloc_sample_array(d, ctx->keyset->bootstrap_key->N);
  res->_signed = _signed;
  res->d = d;
  return res;
}

void ufhe_free_integer(ufhe_integer c){
  free_tlwe_array(c->digits, c->d);
  free(c);
}

void ufhe_copy_integer(ufhe_integer b, ufhe_integer a, ufhe_context ctx){
  for (size_t i = 0; i < MIN(b->d, a->d); i++){
    tlwe_copy(b->digits[i], a->digits[i]);
  }
  if(b->_signed){
    ufhe_extend_integer(b, a->d*ctx->log_torus_base, ctx);
  }
}

void ufhe_encrypt_integer(ufhe_integer out, uint64_t value, ufhe_priv_keyset key, ufhe_context ctx){
  const int mask = ctx->torus_base - 1;
  for (size_t i = 0; i < out->d; i++){
    tlwe_sample(out->digits[i], double2torus(((double)(value&mask))/(ctx->torus_base*2)), key->extracted_trlwe);
    value >>= ctx->log_torus_base;
  }
}

void ufhe_cleartext_integer(ufhe_integer out, uint64_t value, ufhe_context ctx){
  const int mask = ctx->torus_base - 1;
  for (size_t i = 0; i < out->d; i++){
    tlwe_noiseless_trivial_sample(out->digits[i], double2torus(((double)(value&mask))/(ctx->torus_base*2)));
    value >>= ctx->log_torus_base;
  }
}

int64_t ufhe_decrypt_integer(ufhe_integer c, ufhe_priv_keyset key, ufhe_context ctx){
  int64_t result = 0;
  for (int i = c->d - 1; i >= 0; i--){
    result <<= ctx->log_torus_base;
    result |= ((int64_t) round(torus2double(tlwe_phase(c->digits[i], key->extracted_trlwe))*(ctx->torus_base*2))) % ctx->torus_base;
  }
  if(c->_signed) return (result << (64 - ctx->log_torus_base*c->d)) >> (64 - ctx->log_torus_base*c->d);
  else return result;
}

extern ufhe_priv_keyset _priv_key;
extern ufhe_context _ctx;

void ufhe_add_integer(ufhe_integer c, ufhe_integer a, ufhe_integer b, ufhe_context ctx){
  ufhe_sl_add_integer(c, a, 0, b, 0, ctx);
}

void ufhe_extend_integer(ufhe_integer c, int old_precision, ufhe_context ctx){
  const int d_ini = old_precision/ctx->log_torus_base;
  const int n = ctx->keyset->bootstrap_key->n, k = ctx->keyset->bootstrap_key->k, N = ctx->keyset->bootstrap_key->N;
  if(!c->_signed){
    for (size_t i = d_ini; i < c->d; i++) tlwe_noiseless_trivial_sample(c->digits[i], 0);
  }else if(c->d > d_ini){
    TLWE tmp = tlwe_alloc_sample(n);
    TRLWE tmp2 = trlwe_alloc_new_sample(k, N);
    tlwe_keyswitch(tmp, c->digits[d_ini - 1], ctx->keyset->ks_key);
    functional_bootstrap_wo_extract(tmp2, ctx->aux.SIGNEXTEND_LUT, tmp, ctx->keyset->bootstrap_key, ctx->torus_base);
    trlwe_mv_extract_tlwe(&(c->digits[d_ini]), tmp2, c->d - d_ini);
    free_tlwe(tmp);
    free_trlwe(tmp2);
  }
}

/* c = a*B^g + b*B^h */
void ufhe_sl_add_integer(ufhe_integer c, ufhe_integer a, int g, ufhe_integer b, int h, ufhe_context ctx){
  const int n = ctx->keyset->bootstrap_key->n, k = ctx->keyset->bootstrap_key->k, N = ctx->keyset->bootstrap_key->N;
  TLWE tmp = tlwe_alloc_sample(n);
  TRLWE tmp2 = trlwe_alloc_new_sample(k, N);
  const bool _signed = a->_signed || b->_signed;
  const int size = _signed ? a->d : MIN(MAX(a->d + g, b->d + h) + 1, c->d);

  tlwe_noiseless_trivial_sample(c->digits[0], 0);
  for (int i = 0; i < size; i++){
    if(i - g >= 0 && i - g < a->d) tlwe_addto(c->digits[i], a->digits[i - g]);
    if(i - h >= 0 && i - h < b->d) tlwe_addto(c->digits[i], b->digits[i - h]);
    if(i - g < 0 || i - h < 0){ 
      if(i != size - 1) tlwe_noiseless_trivial_sample(c->digits[i + 1], 0);
      continue;
    }
    tlwe_keyswitch(tmp, c->digits[i], ctx->keyset->ks_key);
    functional_bootstrap_wo_extract(tmp2, ctx->aux.ADDSUB_LUT, tmp, ctx->keyset->bootstrap_key, ctx->torus_base);
    trlwe_mv_extract_tlwe_scaling_subto(c->digits[i], tmp2, ctx->torus_base);
    c->digits[i]->b -= double2torus(1./4);
    if(i != size - 1){ // carry
      tlwe_noiseless_trivial_sample(c->digits[i + 1], double2torus(1./(ctx->torus_base*4)));
      trlwe_mv_extract_tlwe_scaling_addto(c->digits[i + 1], tmp2, 1);
    }
  }

  ufhe_extend_integer(c, size*ctx->log_torus_base, ctx);
  free_tlwe(tmp);
  free_trlwe(tmp2);
}

/* b += a*B^g */
void ufhe_sl_addto_integer(ufhe_integer b, ufhe_integer a, int g, ufhe_context ctx){
  const int n = ctx->keyset->bootstrap_key->n, k = ctx->keyset->bootstrap_key->k, N = ctx->keyset->bootstrap_key->N;
  TLWE tmp = tlwe_alloc_sample(n);
  TRLWE tmp2 = trlwe_alloc_new_sample(k, N);
  const bool _signed = a->_signed || b->_signed;
  const int size = _signed ? a->d : MIN(a->d + g + 1, b->d);

  for (int i = 0; i < size; i++){
    if(i - g >= 0 && i - g < a->d) tlwe_addto(b->digits[i], a->digits[i - g]);
    if(i - g < 0) continue;
    tlwe_keyswitch(tmp, b->digits[i], ctx->keyset->ks_key);
    functional_bootstrap_wo_extract(tmp2, ctx->aux.ADDSUB_LUT, tmp, ctx->keyset->bootstrap_key, ctx->torus_base);
    trlwe_mv_extract_tlwe_scaling_subto(b->digits[i], tmp2, ctx->torus_base);
    b->digits[i]->b -= double2torus(1./4);
    if(i != size - 1){ // carry
      b->digits[i + 1]->b += double2torus(1./(ctx->torus_base*4));
      trlwe_mv_extract_tlwe_scaling_addto(b->digits[i + 1], tmp2, 1);
    }
  }

  free_tlwe(tmp);
  free_trlwe(tmp2);
}


void ufhe_sub_integer(ufhe_integer c, ufhe_integer a, ufhe_integer b, ufhe_context ctx){
  const int n = ctx->keyset->bootstrap_key->n, k = ctx->keyset->bootstrap_key->k, N = ctx->keyset->bootstrap_key->N;
  TLWE tmp = tlwe_alloc_sample(n);
  TRLWE tmp2 = trlwe_alloc_new_sample(k, N);

  tlwe_noiseless_trivial_sample(c->digits[0], 0);
  for (int i = 0; i < c->d; i++){
    if(i < a->d) tlwe_addto(c->digits[i], a->digits[i]);
    if(i < b->d) tlwe_subto(c->digits[i], b->digits[i]);
    tlwe_keyswitch(tmp, c->digits[i], ctx->keyset->ks_key);
    functional_bootstrap_wo_extract(tmp2, ctx->aux.ADDSUB_LUT, tmp, ctx->keyset->bootstrap_key, ctx->torus_base);
    trlwe_mv_extract_tlwe_scaling_addto(c->digits[i], tmp2, ctx->torus_base);
    c->digits[i]->b += double2torus(1./4);
    if(i != c->d - 1){ // carry
      tlwe_noiseless_trivial_sample(c->digits[i + 1], double2torus(-1./(ctx->torus_base*4)));
      trlwe_mv_extract_tlwe_scaling_subto(c->digits[i + 1], tmp2, 1);
    }
  }
  free_tlwe(tmp);
  free_trlwe(tmp2);
}

void ufhe_neg_integer(ufhe_integer b, ufhe_integer a, ufhe_context ctx){
  for (int i = 0; i < b->d; i++){ // TODO: sign extension
    tlwe_negate(b->digits[i], a->digits[i]);
    b->digits[i]->b += double2torus(1./2);
    if(i != 0){
      b->digits[i]->b -= double2torus(1./(ctx->torus_base*2));
    }
  }
}

void ufhe_mul_integer(ufhe_integer c, ufhe_integer a, ufhe_integer b, ufhe_context ctx){
  const bool _signed = a->_signed || b->_signed;
  const int size = _signed ? a->d : MIN(a->d + b->d + 1, c->d);

  TLWE sel = ctx->aux.tlwe_selector;
  TLWE * lut_values = ctx->aux.tlwe_lut;

  TRLWE mulmod_lut = ctx->aux.LUT[0];
  TRLWE mulquo_lut = ctx->aux.LUT[1];
  TRLWE * mv_tv = &ctx->aux.LUT[2];

  ufhe_integer prod = ufhe_new_integer(b->d*ctx->log_torus_base, _signed, ctx);
  ufhe_integer carry = ufhe_new_integer(b->d*ctx->log_torus_base, _signed, ctx);
  ufhe_integer res = ufhe_new_integer((b->d + !_signed)*ctx->log_torus_base, _signed, ctx);

  ufhe_cleartext_integer(c, 0, ctx);
  for (size_t i = 0; i < a->d; i++){
    tlwe_keyswitch(sel, a->digits[i], ctx->keyset->ks_key);
    multivalue_bootstrap_phase1(mv_tv, sel, ctx->keyset->bootstrap_key, ctx->torus_base);
    
    // create multiplication mod LUT
    tlwe_noiseless_trivial_sample(lut_values[0], 0);
    tlwe_copy(lut_values[1], a->digits[i]);
    for (size_t j = 2; j < ctx->torus_base; j++){
      multivalue_bootstrap_phase2(lut_values[j], ctx->aux.mulmod_matrix[j], mv_tv, ctx->torus_base, ctx->log_torus_base);
    }
    trlwe_packing_keyswitch(mulmod_lut, lut_values, ctx->keyset->packing_key);
    // create multiplication quo LUT
    tlwe_noiseless_trivial_sample(lut_values[0], 0);
    tlwe_noiseless_trivial_sample(lut_values[1], 0);
    for (size_t j = 2; j < ctx->torus_base; j++){
      multivalue_bootstrap_phase2(lut_values[j], ctx->aux.mulquo_matrix[j], mv_tv, ctx->torus_base, ctx->log_torus_base);
    }
    trlwe_packing_keyswitch(mulquo_lut, lut_values, ctx->keyset->packing_key);
    for (size_t j = 0; j < b->d; j++){
      tlwe_keyswitch(sel, b->digits[j], ctx->keyset->ks_key);
      functional_bootstrap(prod->digits[j], mulmod_lut, sel, ctx->keyset->bootstrap_key, ctx->torus_base);
      functional_bootstrap(carry->digits[j], mulquo_lut, sel, ctx->keyset->bootstrap_key, ctx->torus_base);
      if(i + j >= size) break;
    }
    ufhe_sl_add_integer(res, prod, 0, carry, 1, ctx);
    ufhe_sl_addto_integer(c, res, i, ctx);
  }

  if(c->_signed) ufhe_extend_integer(c, size*ctx->log_torus_base, ctx);
  ufhe_free_integer(prod);
  ufhe_free_integer(carry);
  ufhe_free_integer(res);
}

void ufhe_cmp_integer(ufhe_integer c, ufhe_integer a, ufhe_integer b, ufhe_context ctx){
  const int n = ctx->keyset->bootstrap_key->n, N = ctx->keyset->bootstrap_key->N;
  const int size = MAX(a->d, b->d);

  TLWE tmp_n = tlwe_alloc_sample(n);
  TLWE tmp_N = tlwe_alloc_sample(N);
  TLWE one = tlwe_new_noiseless_trivial_sample(double2torus(1./(2*ctx->torus_base)), N);
  TLWE * lut = (TLWE *) safe_malloc(sizeof(TLWE) * ctx->torus_base);

  ufhe_cleartext_integer(c, 0, ctx);
  lut[0] = c->digits[0];
  for (size_t i = 1; i < ctx->torus_base; i++) lut[i] = one;

  for (size_t i = 0; i < size; i++){
    if(i < a->d && i < b->d){
      tlwe_sub(tmp_N, a->digits[i], b->digits[i]);
      tlwe_keyswitch(tmp_n, tmp_N, ctx->keyset->ks_key);
    }else if(i < a->d){
      tlwe_keyswitch(tmp_n, a->digits[i], ctx->keyset->ks_key);
    }else{
      tlwe_negate(tmp_N, b->digits[i]);
      tlwe_keyswitch(tmp_n, tmp_N, ctx->keyset->ks_key);
    }
    trlwe_packing_keyswitch(ctx->aux.COMPARE_LUT, lut, ctx->keyset->packing_key);
    functional_bootstrap(c->digits[0], ctx->aux.COMPARE_LUT, tmp_n, ctx->keyset->bootstrap_key, ctx->torus_base);
  }
  // signed
  for (size_t i = 0; i < ctx->torus_base/2; i++) lut[i] = c->digits[0];

  if(a->_signed){
    tlwe_negate(one, c->digits[0]);
    trlwe_packing_keyswitch(ctx->aux.COMPARE_LUT, lut, ctx->keyset->packing_key);
    tlwe_keyswitch(tmp_n, a->digits[a->d - 1], ctx->keyset->ks_key);
    functional_bootstrap(c->digits[0], ctx->aux.COMPARE_LUT, tmp_n, ctx->keyset->bootstrap_key, ctx->torus_base);
  }

  if(b->_signed){
    tlwe_negate(one, c->digits[0]);
    trlwe_packing_keyswitch(ctx->aux.COMPARE_LUT, lut, ctx->keyset->packing_key);
    tlwe_keyswitch(tmp_n, b->digits[b->d - 1], ctx->keyset->ks_key);
    functional_bootstrap(c->digits[0], ctx->aux.COMPARE_LUT, tmp_n, ctx->keyset->bootstrap_key, ctx->torus_base);
  }

  c->digits[0]->b += double2torus(1./(2*ctx->torus_base));
  free_tlwe(tmp_n);
  free_tlwe(tmp_N);
  free_tlwe(one);
  free(lut);
}
