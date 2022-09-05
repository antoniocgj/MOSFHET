// This code evaluates a lookup table in leveled mode, encrypting the input as an RGSW sample. 
// A fully homomorphic version requires circuit bootstraps to obtain RGSW samples from (R)LWE samples.
// Notice that some error is expected, as we are trying to get 10-bit precision from very small parameters. Increase the parameters or decrease the output_precision to avoid errors.

#include <mosfhet.h>

TRGSW_DFT encrypt(int m, TRGSW_Key k){
  const uint32_t N = k->trlwe_key->s[0]->N;
  // temporary
  static TRGSW encrypted_input = NULL;
  if(!encrypted_input) encrypted_input = trgsw_alloc_new_sample(k->l, k->Bg_bit, k->trlwe_key->k, N); 
  // encrypting X^(N - m)
  trgsw_monomial_sample(encrypted_input, 1, 2*N - m, k); 
  TRGSW_DFT res = trgsw_alloc_new_DFT_sample(k->l, k->Bg_bit, k->trlwe_key->k, N);
  trgsw_to_DFT(res, encrypted_input);
  return res;
}

void eval_LUT(TLWE output, TRGSW_DFT input, TRLWE LUT){
  // allocating temporaries
  static TRLWE_DFT out_dft = NULL;
  static TRLWE out = NULL;
  if(!out_dft) out_dft = trlwe_alloc_new_DFT_sample(LUT->k, LUT->b->N);
  if(!out) out = trlwe_alloc_new_sample(LUT->k, LUT->b->N);
  // evaluate LUT
  trgsw_mul_trlwe_DFT(out_dft, LUT, input);
  trlwe_from_DFT(out, out_dft);
  trlwe_extract_tlwe(output, out, 0);
}

uint32_t gen_random_input(uint32_t N){
  uint32_t rnd;
  generate_random_bytes(4, (uint8_t *)&rnd);
  return rnd&(N-1); // return random mod N
}

uint32_t * gen_random_LUT(uint32_t size, uint32_t precision){
  uint32_t * res = (uint32_t *) safe_malloc(sizeof(uint32_t)*size);
  generate_random_bytes(sizeof(uint32_t)*size, (uint8_t *)res);
  for (size_t i = 0; i < size; i++){
    res[i] &= (1 << precision) - 1;
  }
  return res;
}

int main(int argc, char const *argv[])
{
  // parameters (for around 7 bits of output precision)
  const uint32_t N = 1024, k = 1, Bg_bit = 8, l = 2;
  const double rlwe_std_dev = 2.989040792967434E-8;
  // setup keys
  TRLWE_Key rlwe_key = trlwe_new_binary_key(N, k, rlwe_std_dev);
  TRGSW_Key input_key = trgsw_new_key(rlwe_key, l, Bg_bit);
  TLWE_Key output_key = tlwe_alloc_key(N, rlwe_std_dev);
  trlwe_extract_tlwe_key(output_key, input_key->trlwe_key);
  // encrypt input
  uint32_t input = gen_random_input(N);
  TRGSW_DFT encrypted_input = encrypt(input, input_key);
  // encrypt lut
  const int output_precision = 10;
  uint32_t * LUT = gen_random_LUT(N, output_precision);
  TRLWE encrypted_LUT = trlwe_new_sample(NULL, rlwe_key);
  for (size_t i = 0; i < N; i++){
    encrypted_LUT->b->coeffs[i] += int2torus(LUT[i], output_precision);
  }
  // eval LUT
  TLWE encrypted_result = tlwe_alloc_sample(N);
  eval_LUT(encrypted_result, encrypted_input, encrypted_LUT);
  // decrypt
  uint32_t result = torus2int(tlwe_phase(encrypted_result, output_key), output_precision);
  // Print
  printf("Expected: %u Result: %u \n", LUT[input], result);
  // frees
  free_trlwe_key(rlwe_key);
  free_trgsw_key(input_key);
  free_tlwe_key(output_key);
  free_trgsw(encrypted_input);
  free_trlwe(encrypted_LUT);
  free_tlwe(encrypted_result);
  // temporary values used in 'eval_LUT' and 'encrypt' should also be freed
  return 0;
}
