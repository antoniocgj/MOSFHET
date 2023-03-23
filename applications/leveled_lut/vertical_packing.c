// CGGI20 Vertical packing LUT evaluation. 
// This code evaluates a big lookup table in leveled mode, encrypting the input bit-by-bit, each one as an RGSW sample. 
// A fully homomorphic version requires circuit bootstraps to obtain RGSW samples from (R)LWE samples.
// Increase the parameters or decrease the output_precision to avoid errors. (or evaluate multiple LUTs)

#include <mosfhet.h>

TRGSW_DFT * encrypt_bits(int m, int size, TRGSW_Key k){
  const uint32_t N = k->trlwe_key->s[0]->N;
  // temporary
  static TRGSW encrypted_input = NULL;
  if(!encrypted_input) encrypted_input = trgsw_alloc_new_sample(k->l, k->Bg_bit, k->trlwe_key->k, N); 
  // encrypting m bit by bit
  TRGSW_DFT * res = trgsw_alloc_new_DFT_sample_array(size, k->l, k->Bg_bit, k->trlwe_key->k, N);
  // m = 2*N - m;
  for (size_t i = 0; i < size; i++){
    trgsw_monomial_sample(encrypted_input, m&1, 0, k); 
    trgsw_to_DFT(res[i], encrypted_input);
    m >>= 1;    
  }
  return res;
}

void CMUX(TRLWE out, TRLWE in1, TRLWE in2, TRGSW_DFT selector){
  TRLWE_DFT tmp = trlwe_alloc_new_DFT_sample(out->k, out->b->N);
  TRLWE tmp2 = trlwe_alloc_new_sample(out->k, out->b->N);
  trlwe_sub(tmp2, in2, in1);
  trgsw_mul_trlwe_DFT(tmp, tmp2, selector);
  trlwe_from_DFT(tmp2, tmp);
  trlwe_add(out, tmp2, in1);
  free_trlwe(tmp);
  free_trlwe(tmp2);
}

// notice that this function destroys LUT
void eval_LUT(TLWE output, TRGSW_DFT * input, int size, TRLWE * LUT){
  const uint32_t N = LUT[0]->b->N, N2 = 2*N, log_N2 = (int) log2(2*N), log_N = log_N2 - 1;  
  // vertical packing
  for (size_t i = 0; i < size - log_N; i++){
    const uint64_t N_luts_half = 1<<(size - log_N - i - 1);
    for (size_t j = 0; j < N_luts_half; j++){
      CMUX(LUT[j], LUT[j], LUT[j + N_luts_half], input[size - i - 1]);
    }
  }
  // evaluate last log_2(N) bits using a blind rotate.
  if(size > log_N) size = log_N;
  // create array of powers of two
  Torus a[32];
  for (size_t i = 0; i < size; i++) a[i] = int2torus(N2 - (1<<i), log_N2);
  blind_rotate(LUT[0], a, input, size);
  trlwe_extract_tlwe(output, LUT[0], 0);
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
  // parameters (for around 5 bits of output precision)
  const uint32_t N = 2048, k = 1, Bg_bit = 23, l = 1;
  const double rlwe_std_dev = 2.2148688116005568e-16;
  const int input_precision = 20, output_precision = 16;
  // setup keys
  TRLWE_Key rlwe_key = trlwe_new_binary_key(N, k, rlwe_std_dev);
  TRGSW_Key input_key = trgsw_new_key(rlwe_key, l, Bg_bit);
  TLWE_Key output_key = tlwe_alloc_key(N, rlwe_std_dev);
  trlwe_extract_tlwe_key(output_key, input_key->trlwe_key);
  // encrypt input
  uint32_t input = gen_random_input(N);
  TRGSW_DFT * encrypted_input = encrypt_bits(input, input_precision, input_key);
  // encrypt lut
  uint32_t * LUT = gen_random_LUT(1<<input_precision, output_precision);
  const uint64_t n_luts = (1<<input_precision)/N;
  TRLWE * encrypted_LUT = trlwe_alloc_new_sample_array(n_luts, k, N);
  for (size_t i = 0; i < n_luts; i++){
    trlwe_sample(encrypted_LUT[i], NULL, rlwe_key);
    for (size_t j = 0; j < N; j++){
      encrypted_LUT[i]->b->coeffs[j] += int2torus(LUT[i*N + j], output_precision);
    }
  }
  // eval LUT
  TLWE encrypted_result = tlwe_alloc_sample(N);
  eval_LUT(encrypted_result, encrypted_input, input_precision, encrypted_LUT);
  // decrypt
  uint32_t result = torus2int(tlwe_phase(encrypted_result, output_key), output_precision);
  // Print
  printf("Expected: %u Result: %u \n", LUT[input], result);
  // frees
  free_trlwe_key(rlwe_key);
  free_trgsw_key(input_key);
  free_tlwe_key(output_key);
  free_trgsw_array(encrypted_input, input_precision);
  free_trlwe_array(encrypted_LUT, n_luts);
  free_tlwe(encrypted_result);
  // temporary values used in 'eval_LUT' and 'encrypt' should also be freed
  return 0;
}
