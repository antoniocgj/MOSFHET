#include "mosfhet.h"
// Functions for bootstrapping using Galois Automorphisms


Bootstrap_GA_Key new_bootstrap_key_ga(TRGSW_Key out_key, TLWE_Key in_key){
  const int l = out_key->l, Bg_bit = out_key->Bg_bit, k = out_key->trlwe_key->k, N = out_key->trlwe_key->s[0]->N;
  Bootstrap_GA_Key res;
  res = (Bootstrap_GA_Key) safe_malloc(sizeof(*res));
  res->s = trgsw_alloc_new_DFT_sample_array(in_key->n, l, Bg_bit, k, N);
  res->ak = trlwe_new_automorphism_KS_keyset(out_key->trlwe_key, true, out_key->l, out_key->Bg_bit);
  res->n = in_key->n;
  res->k = k;
  res->l = l;
  res->N = N;
  res->Bg_bit = Bg_bit;
  res->unfolding = 1;
  TRGSW tmp = trgsw_alloc_new_sample(l, Bg_bit, k, N);
  for (size_t i = 0; i < in_key->n; i++){
    trgsw_monomial_sample(tmp, 1, in_key->s[i], out_key);
    trgsw_to_DFT(res->s[i], tmp);
  }
  free_trgsw(tmp);
  return res;
}

void free_bootstrap_key_ga(Bootstrap_GA_Key key){
  for (size_t i = 0; i < key->n; i++) free_trgsw(key->s[i]);
  free(key->s);
  for (size_t i = 0; i < key->N; i++) free_trlwe_ks_key(key->ak[i]);
  free(key->ak);
  free(key);
}



// Blind rotate 
// Algorithm 4 in https://eprint.iacr.org/2022/198.pdf
// Forcing the all-odd case (since the conversion Torus to integer is approximate anyway)
void blind_rotate_ga(TRLWE tv, Torus * a, TRGSW_DFT * s, TRLWE_KS_Key * ak, int size){
  const uint64_t N = tv->b->N, log_N2 = (uint64_t) log2(2*N), mod_mask = (N<<1) - 1;
  TRLWE rotated_tv = trlwe_alloc_new_sample(tv->k, N);
  TRLWE_DFT tmp = trlwe_alloc_new_DFT_sample(tv->k, N);

  const uint64_t w0p = inverse_mod_2N(torus2int(a[0], log_N2)|1, N);
  trlwe_eval_automorphism(rotated_tv, tv, w0p, ak[(w0p - 1)>>1]);
  for (size_t i = 0; i < size - 1; i++){
    const uint64_t a_i = torus2int(a[i], log_N2)|1;
    const uint64_t w_ip1 = inverse_mod_2N(torus2int(a[i+1], log_N2)|1, N);
    const uint64_t gen = (a_i*w_ip1)&mod_mask;
    trgsw_mul_trlwe_DFT(tmp, rotated_tv, s[i]);
    trlwe_from_DFT(tv, tmp);
    trlwe_eval_automorphism(rotated_tv, tv, gen, ak[(gen - 1)>>1]);
  }
  const uint64_t a_n = torus2int(a[size -1], log_N2)|1;
  trgsw_mul_trlwe_DFT(tmp, rotated_tv, s[size - 1]);
  trlwe_from_DFT(rotated_tv, tmp);
  trlwe_eval_automorphism(tv, rotated_tv, a_n, ak[(a_n - 1)>>1]);
  free_trlwe(rotated_tv);
  free_trlwe(tmp);
}

void functional_bootstrap_wo_extract_ga(TRLWE out, TRLWE tv, TLWE in, Bootstrap_GA_Key key, int torus_base){
  const int N = tv->b->N, N2 = N*2, log_N2 = (int) log2(N*2);
  const Torus prec_offset = double2torus(1./(4*torus_base));
  trlwe_mul_by_xai(out, tv, N2 - torus2int(in->b + prec_offset, log_N2));
  if(key->unfolding == 1) blind_rotate_ga(out, in->a, key->s, key->ak, key->n);
  // else blind_rotate_unfolded(out, in->a, key->su, in->n, key->unfolding); // unimplemented
}

void functional_bootstrap_ga(TLWE out, TRLWE tv, TLWE in, Bootstrap_GA_Key key, int torus_base){
  const int N = tv->b->N;
  TRLWE rotated_tv = trlwe_alloc_new_sample(tv->k, N);
  functional_bootstrap_wo_extract_ga(rotated_tv, tv, in, key, torus_base);
  trlwe_extract_tlwe(out, rotated_tv, 0);
  free_trlwe(rotated_tv);
}