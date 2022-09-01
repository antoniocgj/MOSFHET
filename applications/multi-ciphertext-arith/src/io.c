#include "ufhe.h"

ufhe_priv_keyset ufhe_load_priv_keyset(FILE * fd){
  ufhe_priv_keyset res;
  res = (ufhe_priv_keyset) safe_malloc(sizeof(*res));
  res->tlwe = tlwe_load_new_key(fd);
  res->trlwe = trlwe_load_new_key(fd);
  res->extracted_trlwe = tlwe_alloc_key(res->trlwe->k*res->trlwe->s[0]->N, res->trlwe->sigma);
  trlwe_extract_tlwe_key(res->extracted_trlwe, res->trlwe);
  res->trgsw = trgsw_load_new_key(fd);
  return res;
}

void ufhe_save_priv_keyset(FILE * fd, ufhe_priv_keyset keyset){
  tlwe_save_key(fd, keyset->tlwe);
  trlwe_save_key(fd, keyset->trlwe);
  trgsw_save_key(fd, keyset->trgsw);
}

ufhe_public_keyset ufhe_load_public_keyset(FILE * fd){
  ufhe_public_keyset res;
  res = (ufhe_public_keyset) safe_malloc(sizeof(*res));
  res->bootstrap_key = load_new_bootstrap_key(fd);
  res->ks_key = tlwe_load_new_KS_key(fd);
  res->packing_key = trlwe_load_new_packing_KS_key(fd);
  return res;
}

void ufhe_save_public_keyset(FILE * fd, ufhe_public_keyset keyset){
  save_bootstrap_key(fd, keyset->bootstrap_key);
  tlwe_save_KS_key(fd, keyset->ks_key);
  trlwe_save_packing_KS_key(fd, keyset->packing_key);
}