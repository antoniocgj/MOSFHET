#include "ufhe.h"

// TODO: 2d-table lookup

/* lut[0] = lut[selector] */
void ufhe_encrypted_tlwe_lut(TLWE * selector, TLWE * lut, int size, ufhe_context ctx){
  const int n = ctx->keyset->bootstrap_key->n, k = ctx->keyset->bootstrap_key->k, N = ctx->keyset->bootstrap_key->N;
  TRLWE lut_poly = trlwe_alloc_new_sample(k, N);
  TLWE tmp = tlwe_alloc_sample(n);
  for (size_t i = 0; size > 1; i++){
    tlwe_keyswitch(tmp, selector[i], ctx->keyset->ks_key);
    for (size_t j = 0; j < size/ctx->torus_base; j++){
      trlwe_packing_keyswitch(lut_poly, &lut[j*ctx->torus_base], ctx->keyset->packing_key);
      functional_bootstrap(lut[j], lut_poly, tmp, ctx->keyset->bootstrap_key, ctx->torus_base);
    }
    size /= ctx->torus_base;
  }
  free_trlwe(lut_poly);
  free_tlwe(tmp);
}

/* Evaluate cleartext integer LUT */ 
void ufhe_lut_integer(ufhe_integer out, ufhe_integer selector, uint64_t * lut, int size, ufhe_context ctx){
  const int n = ctx->keyset->bootstrap_key->n, k = ctx->keyset->bootstrap_key->k, N = ctx->keyset->bootstrap_key->N;
  const int mask = ctx->torus_base - 1;
  TLWE * enc_lut = tlwe_alloc_sample_array(size/ctx->torus_base, N);
  TLWE tmp = tlwe_alloc_sample(n);
  TRLWE * c_0_lu = trlwe_alloc_new_sample_array(ctx->torus_base + 1, k, N);
  int * dec_lut = (int *) safe_malloc(sizeof(int) * ctx->torus_base);
  tlwe_keyswitch(tmp, selector->digits[0], ctx->keyset->ks_key);
  multivalue_bootstrap_phase1(c_0_lu, tmp, ctx->keyset->bootstrap_key, ctx->torus_base);
  for (size_t j = 0; j < out->d; j++){
    for (size_t i = 0; i < size/ctx->torus_base; i++){
      for (size_t q = 0; q < ctx->torus_base; q++){
        dec_lut[q] = (lut[i*ctx->torus_base + q] >> (j*ctx->log_torus_base))&mask;
      }
      multivalue_bootstrap_phase2(enc_lut[i], dec_lut, c_0_lu, ctx->torus_base, ctx->log_torus_base); 
    }
    ufhe_encrypted_tlwe_lut(&selector->digits[1], enc_lut, size/ctx->torus_base, ctx);
    tlwe_copy(out->digits[j], enc_lut[0]);
  }
  free_tlwe(tmp);
  free_trlwe_array(c_0_lu, ctx->torus_base +1);
  free(dec_lut);
  free_tlwe_array(enc_lut, size/ctx->torus_base);
}

void ufhe_mux_integer_array(ufhe_integer out, ufhe_integer selector, ufhe_context ctx, int size, ufhe_integer * vec){
  const int N = ctx->keyset->bootstrap_key->N;
  TLWE * enc_lut = tlwe_alloc_sample_array(size, N);
  for (size_t i = 0; i < out->d; i++){
    for (size_t j = 0; j < size; j++){
      tlwe_copy(enc_lut[j], vec[j]->digits[i]);
    }
    if(size%ctx->torus_base){
      for (size_t j = size; j < ctx->torus_base * (1 + size/ctx->torus_base); j++){
        tlwe_copy(enc_lut[j], ctx->aux.tlwe_zero);
      }
    }
    ufhe_encrypted_tlwe_lut(selector->digits, enc_lut, size, ctx);
    tlwe_copy(out->digits[i], enc_lut[0]);
  }
  free_tlwe_array(enc_lut, size);
}
