#include "ufhe.h"

typedef struct{
  // LWE params
  const int n;
  const double lwe_std_dev;
  // RLWE params
  const int N, k;
  const double rlwe_std_dev;
  // RGSW params
  const int l, Bg_bit ;
  // KS params
  const int t, base_bit;
  // PBS params
  const int torus_base;
} mosfhet_params;

static mosfhet_params mosfhet_params_list[3] = {{630, 3.0517578125e-05, 2048, 1, 5.684341886080802e-14, 6, 7, 6, 2, 4}, 
                                            {630, 3.0517578125e-05, 2048, 1, 5.684341886080802e-14, 6, 7, 6, 2, 4}, 
                                            {630, 3.0517578125e-05, 2048, 1, 5.684341886080802e-14, 6, 7, 6, 2, 4}};

ufhe_priv_keyset ufhe_new_priv_keyset(ufhe_params params){
  const mosfhet_params prms = mosfhet_params_list[params];
  ufhe_priv_keyset res;
  res = (ufhe_priv_keyset) safe_malloc(sizeof(*res));
  res->tlwe = tlwe_new_binary_key(prms.n, prms.lwe_std_dev);
  res->trlwe = trlwe_new_binary_key(prms.N, prms.k, prms.rlwe_std_dev);
  res->extracted_trlwe = tlwe_new_binary_key(prms.k*prms.N, prms.lwe_std_dev);
  trlwe_extract_tlwe_key(res->extracted_trlwe, res->trlwe);
  res->trgsw = trgsw_new_key(res->trlwe, prms.l, prms.Bg_bit);
  return res;
}

ufhe_public_keyset ufhe_new_public_keyset(ufhe_priv_keyset priv_key, ufhe_params params){
  const mosfhet_params prms = mosfhet_params_list[params];
  ufhe_public_keyset res;
  res = (ufhe_public_keyset) safe_malloc(sizeof(*res));
  res->bootstrap_key = new_bootstrap_key(priv_key->trgsw, priv_key->tlwe, 1);
  res->ks_key = tlwe_new_KS_key(priv_key->tlwe, priv_key->extracted_trlwe, prms.t, prms.base_bit);
  res->packing_key = trlwe_new_packing_KS_key(priv_key->trlwe, priv_key->extracted_trlwe, prms.t, prms.base_bit, prms.torus_base);
  return res;
}

ufhe_context ufhe_setup_context(ufhe_public_keyset keyset){
  ufhe_context res;
  res = (ufhe_context) safe_malloc(sizeof(*res));
  res->keyset = keyset;
  res->torus_base = keyset->packing_key->torus_base;
  res->log_torus_base = (int) log2(res->torus_base);

  // create constants
  res->aux.tlwe_zero = tlwe_new_noiseless_trivial_sample(0, keyset->bootstrap_key->N);

  // Create LUT slots
  res->aux.lut_torus_slots = (Torus *) safe_malloc(sizeof(Torus)*res->torus_base);
  res->aux.lut_tlwe_slots = (TLWE *) safe_malloc(sizeof(TLWE)*res->torus_base);
  // ADD/SUB LUT [-1/2, -1/2, ..., -1/2]
  res->aux.ADDSUB_LUT = trlwe_alloc_new_sample(keyset->bootstrap_key->k, keyset->bootstrap_key->N);
  Torus tmp = double2torus(-1./(4*res->torus_base));
  trlwe_torus_packing(res->aux.ADDSUB_LUT, &tmp, 1);

  // SIGNEXTEND_LUT [0, 0, ..., 0, B-1, B-1, ..., B-1]
  res->aux.SIGNEXTEND_LUT = trlwe_alloc_new_sample(keyset->bootstrap_key->k, keyset->bootstrap_key->N);
  for (size_t i = 0; i < res->torus_base/2; i++) res->aux.lut_torus_slots[i] = 0;
  for (size_t i = res->torus_base/2; i < res->torus_base; i++) 
    res->aux.lut_torus_slots[i] = double2torus(((double)(res->torus_base - 1))/(2*res->torus_base));
  trlwe_torus_packing(res->aux.SIGNEXTEND_LUT, res->aux.lut_torus_slots, res->torus_base);

  // Compare LUT [0, 1, 1, ..., 1] doesn't require initialization
  res->aux.COMPARE_LUT = trlwe_alloc_new_sample(keyset->bootstrap_key->k, keyset->bootstrap_key->N);

  // Generic 
  #define NUM_OF_TMPS 10
  res->aux.LUT = (TRLWE *) safe_malloc(sizeof(TRLWE) * NUM_OF_TMPS);
  for (size_t i = 0; i < NUM_OF_TMPS; i++){
    res->aux.LUT[i] = trlwe_alloc_new_sample(keyset->bootstrap_key->k, keyset->bootstrap_key->N);
  }
  res->aux.tlwe_selector = tlwe_alloc_sample(keyset->bootstrap_key->n);
  res->aux.tlwe_lut = tlwe_alloc_sample_array(res->torus_base, keyset->bootstrap_key->N);

  // Mul Matrix
  res->aux.mulmod_matrix = (int **) safe_malloc(sizeof(int *)*res->torus_base);
  res->aux.mulquo_matrix = (int **) safe_malloc(sizeof(int *)*res->torus_base);
  for (size_t i = 0; i < res->torus_base; i++){
    res->aux.mulmod_matrix[i] = (int *) safe_malloc(sizeof(int)*res->torus_base);
    res->aux.mulquo_matrix[i] = (int *) safe_malloc(sizeof(int)*res->torus_base);
    for (size_t j = 0; j < res->torus_base; j++){
      res->aux.mulmod_matrix[i][j] = (i*j)%res->torus_base;
      res->aux.mulquo_matrix[i][j] = (i*j)/res->torus_base;
    }
  }
  
  return res;
}