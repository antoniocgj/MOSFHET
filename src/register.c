#include "mosfhet.h"

void free_trgsw_reg(TRGSW_REG p){
  free_trgsw(p->negative);
  free_trgsw(p->positive);
  free(p);
}

void free_trgsw_reg_array(TRGSW_REG * p, int count){
  for (size_t i = 0; i < count; i++){
    free_trgsw_reg(p[i]);
  }
  free(p);
}

TRGSW_REG trgsw_reg_alloc(int l, int Bg_bit, int k, int N){
  TRGSW_REG res;
  res = (TRGSW_REG) safe_malloc(sizeof(*res));
  res->positive = trgsw_alloc_new_DFT_sample(l, Bg_bit, k, N);
  res->negative = trgsw_alloc_new_DFT_sample(l, Bg_bit, k, N);
  return res;
}

TRGSW_REG * trgsw_reg_alloc_array(int count, int l, int Bg_bit, int k, int N){
  TRGSW_REG * res;
  res = (TRGSW_REG *) safe_malloc(sizeof(TRGSW_REG)*count);
  for (size_t i = 0; i < count; i++){
    res[i] = trgsw_reg_alloc(l, Bg_bit, k, N);
  }
  return res;
}

// TRGSW_REG reg_new_trgsw_reg_sample(){
//   reg_alloc_trgsw_reg()
// }

void trgsw_reg_sample(TRGSW_REG out, Torus m, TRGSW_Key key){
  const int N = key->trlwe_key->s[0]->N;
  TRGSW tmp = trgsw_new_exp_sample(m, key);
  trgsw_to_DFT(out->positive, tmp);
  trgsw_monomial_sample(tmp, 1, N - m, key);
  trgsw_to_DFT(out->negative, tmp);
  free_trgsw(tmp);
}

void trgsw_reg_add(TRGSW_REG out, TRGSW_REG in1, TRGSW_REG in2){
  trgsw_mul_DFT2(out->positive, in1->positive, in2->positive);
  trgsw_mul_DFT2(out->negative, in1->negative, in2->negative);
}

void trgsw_reg_negate(TRGSW_REG reg){
  TRGSW_DFT tmp;
  tmp = reg->positive;
  reg->positive = reg->negative;
  reg->negative = tmp;
}

void trgsw_reg_copy(TRGSW_REG out, TRGSW_REG in){
  trgsw_DFT_copy(out->positive, in->positive);
  trgsw_DFT_copy(out->negative, in->negative);
}

void trgsw_reg_sub(TRGSW_REG out, TRGSW_REG in1, TRGSW_REG in2){
  trgsw_mul_DFT2(out->positive, in1->positive, in2->negative);
  trgsw_mul_DFT2(out->negative, in1->negative, in2->positive);
}

void trgsw_reg_subto(TRGSW_REG out, TRGSW_REG in){
  trgsw_mul_DFT2(out->positive, out->positive, in->negative);
  trgsw_mul_DFT2(out->negative, out->negative, in->positive);
}

void trgsw_reg_addto(TRGSW_REG out, TRGSW_REG in1){
  trgsw_reg_add(out, out, in1);
}