#include "mosfhet.h"

Bootstrap_Key new_bootstrap_key_wo_unfolding(TRGSW_Key out_key, TLWE_Key in_key){
  const int l = out_key->l, Bg_bit = out_key->Bg_bit, k = out_key->trlwe_key->k, N = out_key->trlwe_key->s[0]->N;
  Bootstrap_Key res;
  res = (Bootstrap_Key) safe_malloc(sizeof(*res));
  res->s = trgsw_alloc_new_DFT_sample_array(in_key->n, l, Bg_bit, k, N);
  res->n = in_key->n;
  res->k = k;
  res->l = l;
  res->N = N;
  res->Bg_bit = Bg_bit;
  res->unfolding = 1;
  TRGSW tmp = trgsw_alloc_new_sample(l, Bg_bit, k, N);
  for (size_t i = 0; i < in_key->n; i++){
    trgsw_monomial_sample(tmp, in_key->s[i], 0, out_key);
    trgsw_to_DFT(res->s[i], tmp);
  }
  free_trgsw(tmp);
  return res;
}

Bootstrap_Key new_bootstrap_key(TRGSW_Key out_key, TLWE_Key in_key, int unfolding){
  if(unfolding == 1) return new_bootstrap_key_wo_unfolding(out_key, in_key);
  const int l = out_key->l, Bg_bit = out_key->Bg_bit, k = out_key->trlwe_key->k, N = out_key->trlwe_key->s[0]->N;
  Bootstrap_Key res;
  res = (Bootstrap_Key) safe_malloc(sizeof(*res));
  res->n = in_key->n;
  res->k = k;
  res->l = l;
  res->N = N;
  res->Bg_bit = Bg_bit;
  res->unfolding = unfolding;
  const int key_exp = 1 << unfolding, final_exp = key_exp / unfolding; // expansion constants
  res->su = trgsw_alloc_new_sample_array(in_key->n*final_exp, l, Bg_bit, k, N);

  for (size_t i = 0; i < in_key->n; i+=unfolding){
    for (size_t j = 0; j < key_exp; j++){
      Binary key = 1;
      for (size_t u = 0, j_ = j; u < unfolding; u++, j_>>=1){
        if(j_&1) key *= in_key->s[i + u];
        else     key *= 1 - in_key->s[i + u];
      }
      trgsw_monomial_sample(res->su[i*final_exp + j], key, 0, out_key);
    }
  }
  return res;
}


void free_bootstrap_key(Bootstrap_Key key){
  if(key->unfolding == 1){
    for (size_t i = 0; i < key->n; i++) free_trgsw(key->s[i]);
    free(key->s);
  }else{
    const int key_size = key->n*(1 << key->unfolding)/key->unfolding;
    for (size_t i = 0; i < key_size; i++) free_trgsw(key->su[i]);
    free(key->su);
  }
  free(key);
}

void save_bootstrap_key(FILE * fd, Bootstrap_Key key){
  fwrite(&key->n, sizeof(int), 1, fd);
  fwrite(&key->l, sizeof(int), 1, fd);
  fwrite(&key->k, sizeof(int), 1, fd);
  fwrite(&key->N, sizeof(int), 1, fd);
  fwrite(&key->Bg_bit, sizeof(int), 1, fd);
  fwrite(&key->unfolding, sizeof(int), 1, fd);
  if(key->unfolding == 1){
    for (size_t i = 0; i < key->n; i++){
      trgsw_save_DFT_sample(fd, key->s[i]);
    }
  }else{
    const int key_size = key->n*(1 << key->unfolding)/key->unfolding;
    for (size_t i = 0; i < key_size; i++){
      trgsw_save_sample(fd, key->su[i]);
    }
  }
}

Bootstrap_Key load_new_bootstrap_key(FILE * fd){
  Bootstrap_Key res;
  res = (Bootstrap_Key) safe_malloc(sizeof(*res));
  fread(&res->n, sizeof(int), 1, fd);
  fread(&res->l, sizeof(int), 1, fd);
  fread(&res->k, sizeof(int), 1, fd);
  fread(&res->N, sizeof(int), 1, fd);
  fread(&res->Bg_bit, sizeof(int), 1, fd);
  fread(&res->unfolding, sizeof(int), 1, fd);
  if(res->unfolding == 1){
    res->s = (TRGSW_DFT *) safe_malloc(sizeof(TRGSW_DFT) * res->n);
    for (size_t i = 0; i < res->n; i++){
      res->s[i] = trgsw_load_new_DFT_sample(fd, res->l, res->Bg_bit, res->k, res->N);
    }
  }else{
    const int key_size = res->n*(1 << res->unfolding)/res->unfolding;
    res->su = (TRGSW *) safe_malloc(sizeof(TRGSW) * key_size);
    for (size_t i = 0; i < key_size; i++){
      res->su[i] = trgsw_load_new_sample(fd, res->l, res->Bg_bit, res->k, res->N);
    }
  }
  return res;
}


void blind_rotate(TRLWE tv, Torus * a, TRGSW_DFT * s, int size){
  const int N = tv->b->N, log_N2 = (int) log2(2*N);
  TRLWE rotated_tv = trlwe_alloc_new_sample(tv->k, N);
  TRLWE_DFT tmp = trlwe_alloc_new_DFT_sample(tv->k, N);

  for (size_t i = 0; i < size; i++){
    const int a_i = torus2int(a[i], log_N2);
    if(!a_i) continue;
    trlwe_mul_by_xai_minus_1(rotated_tv, tv, a_i);
    trgsw_mul_trlwe_DFT(tmp, rotated_tv, s[i]);
    trlwe_from_DFT(rotated_tv, tmp);
    trlwe_addto(tv, rotated_tv);
  }
  free_trlwe(rotated_tv);
  free_trlwe(tmp);
}

void blind_rotate_unfolded(TRLWE tv, Torus * a, TRGSW * s, int size, int unfolding){
  const int N = tv->b->N, log_N2 = (int) log2(2*N);
  TRLWE_DFT tmp = trlwe_alloc_new_DFT_sample(tv->k, N);
  TRGSW xai = trgsw_alloc_new_sample(s[0]->l, s[0]->Bg_bit, s[0]->samples[0]->k, N);
  TRGSW_DFT xai_DFT = trgsw_alloc_new_DFT_sample(s[0]->l, s[0]->Bg_bit, s[0]->samples[0]->k, N);

  const int key_exp = 1 << unfolding, final_exp = key_exp / unfolding; // expansion constants
  for (size_t i = 0; i < size; i+=unfolding){
    trgsw_copy(xai, s[i*final_exp]);
    for (size_t j = 1; j < key_exp; j++){
      int a_i = 0;
      for (size_t u = 0, j_ = j; u < unfolding; u++, j_>>=1){
        if(j_&1) a_i += torus2int(a[i + u], log_N2);
      } 
      if(!a_i) continue;
      trgsw_mul_by_xai_addto(xai, s[i*final_exp + j], a_i);
    }
    trgsw_to_DFT(xai_DFT, xai);
    trgsw_mul_trlwe_DFT(tmp, tv, xai_DFT);
    trlwe_from_DFT(tv, tmp);
  }

  free_trlwe(tmp);
  free_trgsw(xai);
  free_trgsw(xai_DFT);
}

void functional_bootstrap_wo_extract(TRLWE out, TRLWE tv, TLWE in, Bootstrap_Key key, int torus_base){
  const int N = tv->b->N, N2 = N*2, log_N2 = (int) log2(N*2);
  const Torus prec_offset = double2torus(1./(4*torus_base));
  trlwe_mul_by_xai(out, tv, N2 - torus2int(in->b + prec_offset, log_N2));
  if(key->unfolding == 1) blind_rotate(out, in->a, key->s, in->n);
  else blind_rotate_unfolded(out, in->a, key->su, in->n, key->unfolding);
}

void functional_bootstrap(TLWE out, TRLWE tv, TLWE in, Bootstrap_Key key, int torus_base){
  const int N = tv->b->N;
  TRLWE rotated_tv = trlwe_alloc_new_sample(tv->k, N);
  functional_bootstrap_wo_extract(rotated_tv, tv, in, key, torus_base);
  trlwe_extract_tlwe(out, rotated_tv, 0);
  free_trlwe(rotated_tv);
}

void programmable_bootstrap(TLWE out, TRLWE tv, TLWE in, Bootstrap_Key key, int precision, int kappa, int theta){
  const int N = tv->b->N, N2 = N*2, log_N2 = (int) log2(N2), bit_size = sizeof(Torus)*8;
  const Torus rnd_os = 1ULL << (bit_size - log_N2 + theta - 1);
  const Torus theta_mask = ~((1ULL << (bit_size - log_N2 + theta)) - 1);
  TLWE tmp = tlwe_alloc_sample(in->n);
  for (size_t i = 0; i < in->n; i++){
    tmp->a[i] = in->a[i] << kappa;
    tmp->a[i] = (tmp->a[i] + rnd_os) & theta_mask;
  }
  tmp->b = ((in->b << kappa) + rnd_os) & theta_mask;
  functional_bootstrap(out, tv, tmp, key, 1 << (precision - 1));
  free_tlwe(tmp);
}

void multivalue_bootstrap_CLOT21(TLWE * out, TRLWE tv, TLWE in, Bootstrap_Key key, int torus_base, int n_luts){
  const int slot_size = tv->b->N/(n_luts * torus_base);
  TRLWE tmp = trlwe_alloc_new_sample(tv->k, tv->b->N);
  functional_bootstrap_wo_extract(tmp, tv, in, key, torus_base*n_luts);
  for (size_t i = 0; i < n_luts; i++){
    trlwe_extract_tlwe(out[i], tmp, i*slot_size);
  }
  free_trlwe(tmp);
}

void multivalue_bootstrap_phase1(TRLWE * out, TLWE in, Bootstrap_Key key, int torus_base){
  const int N = out[0]->b->N;
  TRLWE tv = trlwe_new_noiseless_trivial_sample(0, out[0]->k, N);
  for (size_t i = 0; i < N; i++) tv->b->coeffs[i] = double2torus(1./(4*torus_base));
  functional_bootstrap_wo_extract(out[0], tv, in, key, torus_base);
  for (size_t i = 1; i < torus_base; i++){
    trlwe_mul_by_xai(out[i], out[0], i*N/torus_base);
  }
  trlwe_mul_by_xai(out[torus_base], out[0], torus_base);
  trlwe_addto(out[torus_base], out[0]);
  free_trlwe(tv);
}

void multivalue_bootstrap_phase2(TLWE out, int * in, TRLWE * rotated_tv, int torus_base, int log_torus_base){
  const int N = rotated_tv[0]->b->N, k = rotated_tv[0]->k;

  tlwe_noiseless_trivial_sample(out, 0);
  TRLWE tmp = trlwe_alloc_new_sample(k, N);
  for (size_t j = 0; j < log_torus_base; j++){
    const int in_over_tv_0 = ((in[0]>>j)&1) + ((in[torus_base - 1]>>j)&1);
    if(in_over_tv_0 == 2) trlwe_copy(tmp, rotated_tv[torus_base]);
    else if(in_over_tv_0 == 1) trlwe_copy(tmp, rotated_tv[0]);
    else if(in_over_tv_0 == -1) trlwe_negate(tmp, rotated_tv[0]);
    else trlwe_noiseless_trivial_sample(tmp, 0);

    for (size_t i = 1; i < torus_base; i++){
      const int in_over_tv_i = ((in[i]>>j)&1) - ((in[i - 1]>>j)&1);
      if(in_over_tv_i == 1) trlwe_addto(tmp, rotated_tv[i]);
      else if(in_over_tv_i == -1) trlwe_subto(tmp, rotated_tv[i]);
    }    
    trlwe_mv_extract_tlwe_scaling_addto(out, tmp, 1<<j);
  }
  free_trlwe(tmp);
}

void blind_rotate_trgsw(TRGSW tv, Torus * a, TRGSW_DFT * s, int size){
  const int l = tv->l, Bg_bit = tv->Bg_bit, N = tv->samples[0]->b->N, k = tv->samples[0]->k, log_N2 = (int) log2(N*2);
  TRGSW rotated_tv = trgsw_alloc_new_sample(l, Bg_bit, k, N);
  TRGSW_DFT tmp = trgsw_alloc_new_DFT_sample(l, Bg_bit, k, N);

  for (size_t i = 0; i < size; i++){
    const int a_i = torus2int(a[i], log_N2);
    if(!a_i) continue;
    trgsw_mul_by_xai_minus_1(rotated_tv, tv, a_i);
    trgsw_mul_DFT(tmp, rotated_tv, s[i]);
    trgsw_from_DFT(rotated_tv, tmp);
    trgsw_addto(tv, rotated_tv);
  }

  free_trgsw(rotated_tv);
  free_trgsw(tmp);
}

void functional_bootstrap_trgsw_phase1(TRGSW_DFT out, TLWE in, Bootstrap_Key key, int torus_base){
  const int N = out->samples[0]->b->N, l = out->l, Bg_bit = out->Bg_bit, k = out->samples[0]->k, N2 = N*2, log_N2 = (int) log2(N*2);
  const Torus prec_offset = double2torus(1./(4*torus_base));
  TRGSW tv = trgsw_new_noiseless_trivial_sample(1, l, Bg_bit, k, N);
  TRGSW tmp = trgsw_alloc_new_sample(l, Bg_bit, k, N);
  trgsw_mul_by_xai(tmp, tv, N2 - torus2int(in->b + prec_offset, log_N2));
  blind_rotate_trgsw(tmp, in->a, key->s, in->n);
  trgsw_to_DFT(out, tmp);
  free_trgsw(tv);
  free_trgsw(tmp);
}

void functional_bootstrap_trgsw_phase2(TLWE out, TRGSW_DFT in, TRLWE tv){
  const int k = tv->k, N = tv->b->N;
  TRLWE_DFT tmp_dft = trlwe_alloc_new_DFT_sample(k, N);
  TRLWE tmp = trlwe_alloc_new_sample(k, N);
  trgsw_mul_trlwe_DFT(tmp_dft, tv, in);
  trlwe_from_DFT(tmp, tmp_dft);
  trlwe_extract_tlwe(out, tmp, 0);
  free_trlwe(tmp_dft);
  free_trlwe(tmp);
}


void circuit_bootstrap(TRGSW out, TLWE in, Bootstrap_Key key, Generic_KS_Key kska, Generic_KS_Key kskb){
  const int bit_len = sizeof(Torus)*8;
  TRLWE tv = trlwe_alloc_new_sample(key->k, key->N);
  TLWE tmp_out = tlwe_alloc_sample(out->samples[0]->b->N);
  for (size_t i = 0; i < out->l; i++){
    Torus _0h[2] = {0, 1UL << (bit_len - (i + 1) * out->Bg_bit)};
    trlwe_torus_packing(tv, _0h, 2);
    functional_bootstrap(tmp_out, tv, in, key, 2);
    trlwe_priv_keyswitch(out->samples[i], tmp_out, kska);
    trlwe_packing1_keyswitch(out->samples[out->l + i], tmp_out, kskb);
  }
  free_trlwe(tv);
  free_tlwe(tmp_out);
}

void circuit_bootstrap_2(TRGSW out, TLWE in, Bootstrap_Key key, Generic_KS_Key kska, Generic_KS_Key kskb){
  const int bit_len = sizeof(Torus)*8, slot_size = key->N/(2*key->l);
  TRLWE tv = trlwe_alloc_new_sample(key->k, key->N);
  TLWE tmp_out = tlwe_alloc_sample(out->samples[0]->b->N);
  TRLWE tmp = trlwe_alloc_new_sample(tv->k, tv->b->N);
  Torus lut[out->l*2];
  for (size_t i = 0; i < out->l; i++){
    lut[i] = 0;
    lut[key->l + i] = 1ULL << (bit_len - (i + 1) * out->Bg_bit);
  }
  trlwe_torus_packing(tv, lut, 2*out->l);
  functional_bootstrap_wo_extract(tmp, tv, in, key, 2*key->l);
  for (size_t i = 0; i < out->l; i++){
    trlwe_extract_tlwe(tmp_out, tmp, i*slot_size);
    trlwe_priv_keyswitch(out->samples[i], tmp_out, kska);
    trlwe_packing1_keyswitch(out->samples[out->l + i], tmp_out, kskb);
  }
  free_trlwe(tv);
  free_trlwe(tmp);
  free_tlwe(tmp_out);
}

/* out = {p0, p1}[selector] */
void public_mux(TRLWE out, TorusPolynomial p0, TorusPolynomial p1, TRLWE_DFT * selector, int l, int Bg_bit){
  TorusPolynomial p = polynomial_new_torus_polynomial(out->b->N);
  TorusPolynomial * p_dec = polynomial_new_array_of_torus_polynomials(out->b->N, l);
  TRLWE_DFT acc = trlwe_alloc_new_DFT_sample(out->k, out->b->N);
  DFT_Polynomial tmp = polynomial_new_DFT_polynomial(out->b->N);
  polynomial_sub_torus_polynomials(p, p1, p0);
  polynomial_decompose(p_dec, p, Bg_bit, l);
  polynomial_torus_to_DFT(tmp, p_dec[0]);
  trlwe_DFT_mul_by_polynomial(acc, selector[0], tmp);
  for (size_t i = 1; i < l; i++){
    polynomial_torus_to_DFT(tmp, p_dec[i]);
    trlwe_DFT_mul_addto_by_polynomial(acc, selector[i], tmp);
  }
  trlwe_from_DFT(out, acc);
  polynomial_addto_torus_polynomial(out->b, p0);
  // free
  free_polynomial(p);
  free_array_of_polynomials((void **) p_dec, l);
  free_trlwe(acc);
  free_polynomial(tmp);
}

void full_domain_functional_bootstrap_KS21(TLWE out, TorusPolynomial tv, TLWE in, Bootstrap_Key key, Generic_KS_Key ksk, int torus_base){
  const int bit_len = sizeof(Torus)*8;
  const int slot_size = key->N/(key->l * torus_base/2);
  TRLWE tmp_trlwe = trlwe_alloc_new_sample(key->k, key->N);
  TRLWE tmp_trlwe2 = trlwe_alloc_new_sample(key->k, key->N);
  TRLWE_DFT * sign_dec = trlwe_alloc_new_DFT_sample_array(key->l, key->k, key->N);
  TLWE tmp = tlwe_alloc_sample(key->N);
  Torus lut[key->l*torus_base/2];
  for (size_t i = 0; i < key->l; i++){
    for (size_t j = 0; j < torus_base/2; j++){
      lut[i*torus_base/2 + j] = -1ULL << (bit_len - (i + 1) * key->Bg_bit - 1);
    }    
  }
  trlwe_torus_packing_many_LUT(tmp_trlwe, lut, torus_base/2, key->l);
  functional_bootstrap_wo_extract(tmp_trlwe2, tmp_trlwe, in, key, key->l*torus_base/2);
  for (size_t i = 0; i < key->l; i++){
    Torus sign = -1ULL << (bit_len - (i + 1) * key->Bg_bit - 1);
    trlwe_extract_tlwe(tmp, tmp_trlwe2, i*slot_size);
    tmp->b -= sign;
    trlwe_packing1_keyswitch(tmp_trlwe, tmp, ksk);
    trlwe_to_DFT(sign_dec[i], tmp_trlwe);
  }

  TorusPolynomial * p = polynomial_new_array_of_torus_polynomials(tv->N/2, 2);
  for (size_t i = 0; i < tv->N/2; i++){
    p[0]->coeffs[i] = tv->coeffs[i];
    p[1]->coeffs[i] = -tv->coeffs[i + tv->N/2];
  }
  public_mux(tmp_trlwe, p[0], p[1], sign_dec, key->l, key->Bg_bit);
  functional_bootstrap(out, tmp_trlwe, in, key, torus_base/2);
  free_trlwe(tmp_trlwe);
  free_trlwe(tmp_trlwe2);
  free_trlwe_array(sign_dec, key->l);
  free_tlwe(tmp);
  free_array_of_polynomials((void *)p, 2);
}

void full_domain_functional_bootstrap_KS21_2(TLWE out, TorusPolynomial tv, TLWE in, Bootstrap_Key key, Generic_KS_Key ksk, int torus_base){
  const int bit_len = sizeof(Torus)*8;
  TRLWE tmp_trlwe = trlwe_alloc_new_sample(key->k, key->N);
  TRLWE_DFT * sign_dec = trlwe_alloc_new_DFT_sample_array(key->l, key->k, key->N);
  TLWE tmp = tlwe_alloc_sample(key->N);
  for (size_t i = 0; i < key->l; i++){
    Torus sign[1] = {-1ULL << (bit_len - (i + 1) * key->Bg_bit - 1)};
    trlwe_torus_packing(tmp_trlwe, sign, 1);
    functional_bootstrap(tmp, tmp_trlwe, in, key, torus_base/2);
    tmp->b -= sign[0];

    trlwe_packing1_keyswitch(tmp_trlwe, tmp, ksk);
    trlwe_to_DFT(sign_dec[i], tmp_trlwe);
  }

  TorusPolynomial * p = polynomial_new_array_of_torus_polynomials(tv->N/2, 2);
  for (size_t i = 0; i < tv->N/2; i++){
    p[0]->coeffs[i] = tv->coeffs[i];
    p[1]->coeffs[i] = -tv->coeffs[i + tv->N/2];
  }
  public_mux(tmp_trlwe, p[0], p[1], sign_dec, key->l, key->Bg_bit);
  functional_bootstrap(out, tmp_trlwe, in, key, torus_base/2);
  free_trlwe(tmp_trlwe);
  free_trlwe_array(sign_dec, key->l);
  free_tlwe(tmp);
  free_array_of_polynomials((void *)p, 2);
}

void full_domain_functional_bootstrap_CLOT21(TLWE out, TRLWE tv[2], TLWE in, Bootstrap_Key key, Generic_KS_Key ksk, TRLWE_KS_Key rlk, int precision){
  const int bit_len = sizeof(Torus)*8;
  TRLWE tmp_trlwe = trlwe_alloc_new_sample(key->k, key->N);
  TLWE ct_sign = tlwe_alloc_sample(key->N);
  TLWE ct_f0 = tlwe_alloc_sample(key->N);
  TLWE ct_f1 = tlwe_alloc_sample(key->N);

  Torus sign[1] = {1ULL << (bit_len - precision - 1)};
  trlwe_torus_packing(tmp_trlwe, sign, 1);

  functional_bootstrap(ct_f0, tv[0], in, key, 1 << (precision - 1));
  functional_bootstrap(ct_f1, tv[1], in, key, 1 << (precision - 1));
  functional_bootstrap(ct_sign, tmp_trlwe, in, key, 1 << (precision - 1));

  ct_sign->b -= sign[0];
  tlwe_mul(ct_f1, ct_f1, ct_sign, precision, ksk, rlk);
  ct_sign->b += 2*sign[0];
  tlwe_mul(ct_f0, ct_f0, ct_sign, precision, ksk, rlk);

  tlwe_add(out, ct_f0, ct_f1);

  free_trlwe(tmp_trlwe);
  free_tlwe(ct_f0);
  free_tlwe(ct_f1);
  free_tlwe(ct_sign);
}

void full_domain_functional_bootstrap_CLOT21_2(TLWE out, Torus * tv, TLWE in, Bootstrap_Key key, Generic_KS_Key ksk, TRLWE_KS_Key rlk, int precision){
  const int bit_len = sizeof(Torus)*8, torus_base = 1 << (precision - 2);
  const int slot_size = key->N/(4 * torus_base);
  TRLWE tmp_trlwe = trlwe_alloc_new_sample(key->k, key->N);
  TRLWE tmp_trlwe2 = trlwe_alloc_new_sample(key->k, key->N);
  TLWE ct_sign = tlwe_alloc_sample(key->N);
  TLWE ct_f0 = tlwe_alloc_sample(key->N);
  TLWE ct_f1 = tlwe_alloc_sample(key->N);

  Torus sign = 1ULL << (bit_len - precision - 1);
  Torus lut[4*torus_base];
  memcpy(lut, tv, sizeof(Torus)*2*torus_base);
  for (size_t i = 2*torus_base; i < 3*torus_base; i++){
    lut[i] = sign;
  }
  trlwe_torus_packing_many_LUT(tmp_trlwe, lut, torus_base, 4);

  functional_bootstrap_wo_extract(tmp_trlwe2, tmp_trlwe, in, key, 4*torus_base);
  trlwe_extract_tlwe(ct_f0, tmp_trlwe2, 0);
  trlwe_extract_tlwe(ct_f1, tmp_trlwe2, slot_size);
  trlwe_extract_tlwe(ct_sign, tmp_trlwe2, 2*slot_size);

  ct_sign->b -= sign;
  tlwe_mul(ct_f1, ct_f1, ct_sign, precision, ksk, rlk);
  ct_sign->b += 2*sign;
  tlwe_mul(ct_f0, ct_f0, ct_sign, precision, ksk, rlk);

  tlwe_add(out, ct_f0, ct_f1);

  free_trlwe(tmp_trlwe);
  free_trlwe(tmp_trlwe2);
  free_tlwe(ct_f0);
  free_tlwe(ct_f1);
  free_tlwe(ct_sign);
}

void full_domain_functional_bootstrap(TLWE out, TRLWE tv, TLWE in, Bootstrap_Key key, TLWE_KS_Key tlwe_ksk, int precision){
  const int bit_len = sizeof(Torus)*8;
  TRLWE tmp_trlwe = trlwe_alloc_new_sample(key->k, key->N);
  TLWE ct_sign = tlwe_alloc_sample(key->N);
  TLWE in2 = tlwe_alloc_sample(in->n);

  Torus sign[1] = {(1ULL << (bit_len - 2)) - (1ULL << (bit_len - precision - 2))};
  trlwe_torus_packing(tmp_trlwe, sign, 1);
  functional_bootstrap(ct_sign, tmp_trlwe, in, key, 1 << (precision - 1));
  ct_sign->b -= sign[0];

  tlwe_keyswitch(in2, ct_sign, tlwe_ksk);
  tlwe_addto(in2, in);
  
  functional_bootstrap(out, tv, in2, key, 1 << (precision));

  free_trlwe(tmp_trlwe);
  free_tlwe(ct_sign);
  free_tlwe(in2);
}

