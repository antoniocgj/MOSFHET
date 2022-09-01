#include "ufhe.h"


void ufhe_relu_integer(ufhe_integer out, ufhe_integer in, ufhe_context ctx){
  tlwe_keyswitch(ctx->aux.tlwe_selector, in->digits[in->d - 1], ctx->keyset->ks_key);
  for (size_t i = 0; i < in->d - 1; i++){
    for (size_t j = 0; j < ctx->torus_base/2; j++){
      ctx->aux.lut_tlwe_slots[j] = in->digits[i];
      ctx->aux.lut_tlwe_slots[j + ctx->torus_base/2] = ctx->aux.tlwe_zero;
    }
    trlwe_packing_keyswitch(ctx->aux.LUT[0], ctx->aux.lut_tlwe_slots, ctx->keyset->packing_key);
    functional_bootstrap(out->digits[i], ctx->aux.LUT[0] , ctx->aux.tlwe_selector, ctx->keyset->bootstrap_key, ctx->torus_base);
  }
  
  for (size_t j = 0; j < ctx->torus_base/2; j++){
    ctx->aux.lut_torus_slots[j] = double2torus(((double) j)/(ctx->torus_base*2));
    ctx->aux.lut_torus_slots[j + ctx->torus_base/2] = double2torus(0);
  }
  trlwe_torus_packing(ctx->aux.LUT[0], ctx->aux.lut_torus_slots, ctx->torus_base);
  functional_bootstrap(out->digits[in->d - 1], ctx->aux.LUT[0], ctx->aux.tlwe_selector, ctx->keyset->bootstrap_key, ctx->torus_base);
}