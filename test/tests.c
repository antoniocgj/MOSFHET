#include "unity_test/unity.h"
#include <mosfhet.h>

// Global debug variables
TRGSW_Key _glb_debug_trgsw_key;
TRLWE_Key _glb_debug_trlwe_key;
TLWE_Key _glb_debug_tlwe_key;
TLWE_Key _glb_debug_tlwe_key2;




// TEST Macros
#ifdef TORUS32
#define TEST_ASSERT_TORUS_WITHIN_MESSAGE(A, B, C, D) TEST_ASSERT_HEX32_WITHIN_MESSAGE(((A)>>32), (B), (C), (D))
#define TEST_ASSERT_TORUS_ARRAY_WITHIN_MESSAGE(A, B, C, D, E) TEST_ASSERT_HEX32_ARRAY_WITHIN_MESSAGE(((A)>>32), (B), (C), (D), (E))
#define SKIP_IF_TORUS32 TEST_IGNORE_MESSAGE("Requires Torus64 or not implemented for Torus32");

// LWE params
const int n = 630;
const double lwe_std_dev = 3.0517578125e-05; // 2^-15
// RLWE params
const int N = 1024, k = 1;
const double rlwe_std_dev = 2.9802322387695312e-08; // 2^-25
// RGSW params
// const int l = 6, Bg_bit = 6;
const int l = 5, Bg_bit = 5;
// KS params
const int t = 2, base_bit = 6;

#else
#define TEST_ASSERT_TORUS_WITHIN_MESSAGE TEST_ASSERT_HEX64_WITHIN_MESSAGE
#define TEST_ASSERT_TORUS_ARRAY_WITHIN_MESSAGE TEST_ASSERT_INT64_ARRAY_WITHIN_MESSAGE
#define SKIP_IF_TORUS32 

// Note: BR Unfolding requires n to be divisible by the unfolding value. Example: In SET_1, you need to change n from 585 to 592, since 8|592.

// From eprint 2022/704 table 4
#define SET_4
// set 1 (set 1 should fail most tests)
#if defined(SET_1)
const int n = 585, N = 1024, k = 1, Bg_bit = 8, l = 2, base_bit = 2, t = 5;
const double lwe_std_dev = 9.141776004202573E-5, rlwe_std_dev = 2.989040792967434E-8;
#elif defined(SET_2)
// set 2
const int n = 720, N = 2048, k = 1, Bg_bit = 23, l = 1, base_bit = 4 , t = 5;
const double lwe_std_dev = 7.747831515176779e-6, rlwe_std_dev = 2.2148688116005568e-16;
#elif defined(SET_3)
// set 3
const int n = 829, N = 4096, k = 1, Bg_bit = 23, l = 1, base_bit = 2, t = 11;
const double lwe_std_dev = 1.0562341599676662e-6, rlwe_std_dev = 2.168404344971009e-19;
#else
// From TFHEpp
// LWE params
const int n = 632;
const double lwe_std_dev = 3.0517578125e-05; // 2^-15
// RLWE params
const int N = 2048, k = 1;
const double rlwe_std_dev = 5.684341886080802e-14; // 2^-44
// RGSW params
// const int l = 6, Bg_bit = 6;
const int l = 4, Bg_bit = 9;
// KS params
const int t = 8, base_bit = 4;
#endif

#endif

void setUp(void) {}
void tearDown(void) {}

void test_normal_generator(){
  double test[100000], var = 0;
  for (size_t i = 0; i < 100000; i++){
    test[i] = generate_normal_random(0.001);
    var += test[i]*test[i];
  }
  var /= 100000 - 1;
  TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(0.00001, 0.000001, var, "Variance of one variable.");
  
  var = 0;
  for (size_t i = 0; i < 100000; i++){
    test[i] = generate_normal_random(0.001) + generate_normal_random(0.001);
    var += test[i]*test[i];
  }
  var /= 100000 - 1;
  TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(0.00001, 0.000002, var, "Variance of adding two variables.");

  var = 0;
  for (size_t i = 0; i < 100000; i++){
    test[i] = 2*generate_normal_random(0.001);
    var += test[i]*test[i];
  }
  var /= 100000 - 1;
  TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(0.00001, 0.000004, var, "Variance of multiplying a variable.");

  var = 0;
  for (size_t i = 0; i < 1000; i++){
    test[i] = 0;
    for (size_t j = 0; j < 10000; j++){
      test[i] += generate_normal_random(0.001);
    }
    var += test[i]*test[i];
  }
  var /= 1000 - 1;
  // printf("var %lf\n", var);
  TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(0.001, 0.01, var, "Variance of adding 10000 variables.");

  var = 0;
  for (size_t i = 0; i < 1000; i++){
    Torus testT = 0;
    for (size_t j = 0; j < 10000; j++){
      const double val = generate_normal_random(0.001);
      testT += double2torus(val);
    }
    if(testT > (1ULL << (sizeof(Torus)*8 - 1))) testT = -1ULL - testT;
    var += torus2double(testT)*torus2double(testT);
  }
  var /= 1000 - 1;
  // printf("var %lf\n", var);
  TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(0.001, 0.01, var, "Variance of adding 10000 Torus variables.");
}

void test_tlwe(){
  TLWE_Key key = tlwe_new_binary_key(n, lwe_std_dev);
  Torus in1, in2;
  generate_random_bytes(sizeof(Torus), (uint8_t *) &in1);
  generate_random_bytes(sizeof(Torus), (uint8_t *) &in2);
  TLWE c_1 = tlwe_new_sample(in1, key);
  TEST_ASSERT_TORUS_WITHIN_MESSAGE(1UL << 54, in1, tlwe_phase(c_1, key), "TLWE decryption failed");
  TLWE c_2 = tlwe_new_sample(in2, key);
  TLWE c_3 = tlwe_new_noiseless_trivial_sample(0, n);
  tlwe_add(c_3, c_1, c_2);
  TEST_ASSERT_TORUS_WITHIN_MESSAGE(1UL << 54, in1+in2, tlwe_phase(c_3, key), "TLWE addition failed");
  tlwe_addto(c_3, c_1);
  TEST_ASSERT_TORUS_WITHIN_MESSAGE(1UL << 54, in1+in2+in1, tlwe_phase(c_3, key), "TLWE addto failed");
  tlwe_sub(c_3, c_1, c_2);
  TEST_ASSERT_TORUS_WITHIN_MESSAGE(1UL << 54, in1-in2, tlwe_phase(c_3, key), "TLWE sub failed");
  tlwe_subto(c_3, c_1);
  TEST_ASSERT_TORUS_WITHIN_MESSAGE(1UL << 54, in1-in2-in1, tlwe_phase(c_3, key), "TLWE subto failed");
  tlwe_noiseless_trivial_sample(c_3, 0);
  in1 = 0;
  for (size_t i = 0; i < 10000; i++){
    generate_random_bytes(sizeof(Torus), (uint8_t *) &in2);
    in1 += in2;
    tlwe_sample(c_2, in2, key);
    tlwe_addto(c_3, c_2);
    // printf("%ld \n", i);
    // const Torus res = tlwe_phase(c_3, key);
    // if(labs(res - in1) > 1UL << 58 && labs(res - in1) < 1UL << 63) break;
  }
  TEST_ASSERT_TORUS_WITHIN_MESSAGE(1UL << 58, in1, tlwe_phase(c_3, key), "10000 TLWE additions failed");

  free_tlwe(c_1);
  free_tlwe(c_2);
  free_tlwe(c_3);
  free_tlwe_key(key);
}

void test_trlwe(){
  TRLWE_Key key = trlwe_new_binary_key(N, k, rlwe_std_dev);
  TorusPolynomial poly_1 = polynomial_new_torus_polynomial(N);
  generate_random_bytes(N*sizeof(Torus), (uint8_t *) poly_1->coeffs);
  
  TorusPolynomial poly_res = polynomial_new_torus_polynomial(N);
  TorusPolynomial poly_sum = polynomial_new_torus_polynomial(N);
  TRLWE c_1 = trlwe_new_sample(poly_1, key);
  trlwe_phase(poly_res, c_1, key);
  TEST_ASSERT_TORUS_ARRAY_WITHIN_MESSAGE(1UL << 44, poly_1->coeffs, poly_res->coeffs, N, "TRLWE decryption failed.");

  TorusPolynomial poly_2 = polynomial_new_torus_polynomial(N);
  generate_random_bytes(N*sizeof(Torus), (uint8_t *) poly_2->coeffs);
  TRLWE c_2 = trlwe_new_sample(poly_2, key);
  TRLWE c_3 = trlwe_new_noiseless_trivial_sample(NULL, k, N);
  trlwe_add(c_3, c_1, c_2);
  polynomial_add_torus_polynomials(poly_sum, poly_1, poly_2);
  trlwe_phase(poly_res, c_3, key);
  TEST_ASSERT_TORUS_ARRAY_WITHIN_MESSAGE(1UL << 44, poly_sum->coeffs, poly_res->coeffs, N, "TRLWE addition failed.");

  trlwe_sub(c_3, c_1, c_2);
  polynomial_sub_torus_polynomials(poly_sum, poly_1, poly_2);
  trlwe_phase(poly_res, c_3, key);
  TEST_ASSERT_TORUS_ARRAY_WITHIN_MESSAGE(1UL << 44, poly_sum->coeffs, poly_res->coeffs, N, "TRLWE subtraction failed.");

  free_polynomial(poly_1);
  free_polynomial(poly_2);
  free_polynomial(poly_sum);
  free_polynomial(poly_res);
  free_trlwe(c_1);
  free_trlwe(c_2);
  free_trlwe(c_3);
  free_trlwe_key(key);
}

void test_compressed_trlwe(){
  SKIP_IF_TORUS32
#ifdef USE_VAES
  aes_setup_rnd_seed();
#endif
  TRLWE_Key key = trlwe_new_binary_key(N, k, rlwe_std_dev);
  TorusPolynomial poly_1 = polynomial_new_torus_polynomial(N);
  generate_random_bytes(N*sizeof(Torus), (uint8_t *) poly_1->coeffs);

  TorusPolynomial poly_2 = polynomial_new_torus_polynomial(N);
  generate_random_bytes(N*sizeof(Torus), (uint8_t *) poly_2->coeffs);
  
  TorusPolynomial poly_3 = polynomial_new_torus_polynomial(N);
  generate_random_bytes(N*sizeof(Torus), (uint8_t *) poly_3->coeffs);
  
  TorusPolynomial poly_res = polynomial_new_torus_polynomial(N);
  TorusPolynomial poly_sum = polynomial_new_torus_polynomial(N);
  TRLWE c_1 = trlwe_new_sample(poly_1, key);
  TRLWE c_2 = trlwe_new_compressed_sample(poly_2, key);
  trlwe_compressed_subto(c_1, c_2);
  polynomial_sub_torus_polynomials(poly_sum, poly_1, poly_2);
  trlwe_phase(poly_res, c_1, key);
  TEST_ASSERT_TORUS_ARRAY_WITHIN_MESSAGE(1UL << 44, poly_sum->coeffs, poly_res->coeffs, N, "Compressed TRLWE subto failed.");

  TRLWE c_3 = trlwe_new_compressed_sample(poly_3, key);
  trlwe_compressed_subto(c_1, c_3);
  polynomial_sub_torus_polynomials(poly_sum, poly_sum, poly_3);
  trlwe_phase(poly_res, c_1, key);
  TEST_ASSERT_TORUS_ARRAY_WITHIN_MESSAGE(1UL << 44, poly_sum->coeffs, poly_res->coeffs, N, "Compressed TRLWE subto failed.");

  free_polynomial(poly_1);
  free_polynomial(poly_2);
  free_polynomial(poly_3);
  free_polynomial(poly_sum);
  free_polynomial(poly_res);
  free_trlwe(c_1);
  free_trlwe(c_2);
  free_trlwe(c_3);
  free_trlwe_key(key);
}

void test_compressed_trlwe_rotate_vaes(){
  SKIP_IF_TORUS32
#ifdef USE_VAES
  aes_setup_rnd_seed();
#endif
  TRLWE_Key key = trlwe_new_binary_key(N, k, rlwe_std_dev);
  TorusPolynomial poly_1 = polynomial_new_torus_polynomial(N);
  generate_random_bytes(N*sizeof(Torus), (uint8_t *) poly_1->coeffs);

  TorusPolynomial poly_2 = polynomial_new_torus_polynomial(N);
  generate_random_bytes(N*sizeof(Torus), (uint8_t *) poly_2->coeffs);
  
  TorusPolynomial poly_3 = polynomial_new_torus_polynomial(N);
  generate_random_bytes(N*sizeof(Torus), (uint8_t *) poly_3->coeffs);
  
  TorusPolynomial poly_res = polynomial_new_torus_polynomial(N);
  TorusPolynomial poly_sum = polynomial_new_torus_polynomial(N);
  TRLWE c_1 = trlwe_new_sample(poly_1, key);
  TRLWE c_2 = trlwe_new_compressed_sample(poly_2, key);
  polynomial_copy_torus_polynomial(poly_sum, c_1->b);
#ifdef USE_VAES
  trlwe_mul_by_xai_addto_comp_vaes(c_1, c_2, 458);
#endif
  torus_polynomial_mul_by_xai_addto(poly_1, poly_2, 458);
  torus_polynomial_mul_by_xai_addto(poly_sum, c_2->b, 458);
  trlwe_phase(poly_res, c_1, key);
  TEST_ASSERT_TORUS_ARRAY_WITHIN_MESSAGE(1UL << 44, poly_sum->coeffs, c_1->b->coeffs, N, "Compressed TRLWE->b mul by xai failed.");
  TEST_ASSERT_TORUS_ARRAY_WITHIN_MESSAGE(1UL << 44, poly_1->coeffs, poly_res->coeffs, N, "Compressed TRLWE mul by xai failed.");

  // TRLWE c_3 = trlwe_new_compressed_sample(poly_3, key);
  // trlwe_compressed_subto(c_1, c_3);
  // polynomial_sub_torus_polynomials(poly_sum, poly_sum, poly_3);
  // trlwe_phase(poly_res, c_1, key);
  // TEST_ASSERT_TORUS_ARRAY_WITHIN_MESSAGE(1UL << 44, poly_sum->coeffs, poly_res->coeffs, N, "Compressed TRLWE subto failed.");

  free_polynomial(poly_1);
  free_polynomial(poly_2);
  free_polynomial(poly_3);
  free_polynomial(poly_sum);
  free_polynomial(poly_res);
  free_trlwe(c_1);
  free_trlwe(c_2);
  // free_trlwe(c_3);
  free_trlwe_key(key);
}

void test_poly_DFT(){
  TorusPolynomial poly = polynomial_new_torus_polynomial(N);
  generate_random_bytes(N*sizeof(Torus), (uint8_t *) poly->coeffs);
  TorusPolynomial poly_res = polynomial_new_torus_polynomial(N);
  DFT_Polynomial poly_dft = polynomial_new_DFT_polynomial(N);
  polynomial_torus_to_DFT(poly_dft, poly);
  polynomial_DFT_to_torus(poly_res, poly_dft);
  TEST_ASSERT_TORUS_ARRAY_WITHIN_MESSAGE(1UL << 40, poly->coeffs, poly_res->coeffs, N, "DFT transform failed.");
  free_polynomial(poly);
  free_polynomial(poly_res);
  free_polynomial(poly_dft);
}

void test_poly_DFT_mul(){
  TorusPolynomial poly_1 = polynomial_new_torus_polynomial(N);
  generate_random_bytes(N*sizeof(Torus), (uint8_t *) poly_1->coeffs);
  TorusPolynomial poly_2 = polynomial_new_torus_polynomial(N);
  generate_random_bytes(N*sizeof(Torus), (uint8_t *) poly_2->coeffs);
  const int scale = sizeof(Torus)*8 - 10;
  for (size_t i = 0; i < N; i++) poly_2->coeffs[i] >>= scale;
  TorusPolynomial poly_res_1 = polynomial_new_torus_polynomial(N);
  TorusPolynomial poly_res_2 = polynomial_new_torus_polynomial(N);
  DFT_Polynomial poly_dft_1 = polynomial_new_DFT_polynomial(N);
  DFT_Polynomial poly_dft_2 = polynomial_new_DFT_polynomial(N);
  DFT_Polynomial poly_dft_3 = polynomial_new_DFT_polynomial(N);
  polynomial_torus_to_DFT(poly_dft_1, poly_1);
  polynomial_torus_to_DFT(poly_dft_2, poly_2);
  polynomial_mul_DFT(poly_dft_3, poly_dft_1, poly_dft_2);
  polynomial_DFT_to_torus(poly_res_1, poly_dft_3);
  memset(poly_res_2->coeffs, 0, N*sizeof(Torus));
  polynomial_naive_mul_addto_torus(poly_res_2, poly_1, poly_2);
  TEST_ASSERT_TORUS_ARRAY_WITHIN_MESSAGE(1UL << 40, poly_res_2->coeffs, poly_res_1->coeffs, N, "DFT multiplication failed.");

  polynomial_mul_addto_DFT(poly_dft_3, poly_dft_1, poly_dft_2);
  polynomial_DFT_to_torus(poly_res_1, poly_dft_3);
  polynomial_naive_mul_addto_torus(poly_res_2, poly_1, poly_2);
  TEST_ASSERT_TORUS_ARRAY_WITHIN_MESSAGE(1UL << 40, poly_res_2->coeffs, poly_res_1->coeffs, N, "DFT mul_addto failed.");
  free_polynomial(poly_1);
  free_polynomial(poly_2);
  free_polynomial(poly_res_1);
  free_polynomial(poly_res_2);
  free_polynomial(poly_dft_1);
  free_polynomial(poly_dft_2);
  free_polynomial(poly_dft_3);
}


static void _debug_print128(char * msg, __uint128_t x){
  uint64_t * x64 = (uint64_t *) &x;
  printf("%s: 0x%016lx%016lx\n", msg, x64[1], x64[0]);
}


// void karatsuba_u128_scale64(uint64_t * restrict out, uint64_t * restrict in1, uint64_t * restrict in2, int size, int bit_scale);
// void poly_mul_int64to128(__uint128_t * restrict out, uint64_t * restrict p, uint64_t * restrict p2, int size);
// void test_poly_int128_mul(){
//   TorusPolynomial poly_1 = polynomial_new_torus_polynomial(N);
//   generate_random_bytes(N*sizeof(Torus), (uint8_t *) poly_1->coeffs);
//   TorusPolynomial poly_2 = polynomial_new_torus_polynomial(N);
//   generate_random_bytes(N*sizeof(Torus), (uint8_t *) poly_2->coeffs);
//   TorusPolynomial poly_3 = polynomial_new_torus_polynomial(N);
//   generate_random_bytes(N*sizeof(Torus), (uint8_t *) poly_3->coeffs);
//   TorusPolynomial poly_4 = polynomial_new_torus_polynomial(N);
//   generate_random_bytes(N*sizeof(Torus), (uint8_t *) poly_4->coeffs);
 
 
//   __uint128_t * tmp = (__uint128_t *)safe_aligned_malloc(sizeof(__uint128_t)*N*2);

//   karatsuba_u128_scale64(poly_3->coeffs, poly_1->coeffs, poly_2->coeffs, N, 60); 
//   poly_mul_int64to128(tmp, poly_1->coeffs, poly_2->coeffs, N);
//   for (size_t i = 0; i < N; i++){
//     poly_4->coeffs[i] = (uint64_t)(tmp[i]>>60) - (uint64_t)(tmp[N + i]>>60);
//   }

//   TEST_ASSERT_TORUS_ARRAY_WITHIN_MESSAGE(1UL << 62, poly_4->coeffs, poly_3->coeffs, N, "Karatsuba u128 mul failed.");

//   free(tmp);
//   free_polynomial(poly_1);
//   free_polynomial(poly_2);
//   free_polynomial(poly_3);
//   free_polynomial(poly_4);
// }


void test_trgsw(){
  TRLWE_Key trlwe_key = trlwe_new_binary_key(N, k, rlwe_std_dev);
  uint64_t test;
  generate_random_bytes(sizeof(uint64_t), (uint8_t *) &test);
  test &= (N - 1);
  TRGSW_Key key = trgsw_new_key(trlwe_key, l, Bg_bit);
  TRGSW c_1 = trgsw_new_exp_sample(test, key);
  uint64_t res = _debug_trgsw_decrypt_exp_sample(c_1, key);
  TEST_ASSERT_EQUAL_INT64_MESSAGE(test, res, "TRGSW EXP ENCRYPTION FAILED.");
  free_trlwe_key(trlwe_key);
  free_trgsw_key(key);
  free_trgsw(c_1);
}


void test_trgsw_sub(){
  TRLWE_Key trlwe_key = trlwe_new_binary_key(N, k, rlwe_std_dev);
  TRGSW_Key key = trgsw_new_key(trlwe_key, l, Bg_bit);

  TRGSW c_1 = trgsw_new_sample(1, key);
  TRGSW c_2 = trgsw_new_sample(0, key);
  TRGSW z_inv = trgsw_alloc_new_sample(l, Bg_bit, k, N), 
        one = trgsw_new_noiseless_trivial_sample(1, l, Bg_bit, k, N);

  trgsw_sub(z_inv, one, c_1);
  uint64_t res = _debug_trgsw_decrypt_exp_sample(z_inv, key);
  TEST_ASSERT_EQUAL_INT64_MESSAGE(-1, res, "TRGSW SUB FAILED.");

  trgsw_sub(z_inv, one, c_2);
  res = _debug_trgsw_decrypt_exp_sample(z_inv, key);
  TEST_ASSERT_EQUAL_INT64_MESSAGE(0, res, "TRGSW SUB FAILED.");
  free_trlwe_key(trlwe_key);
  free_trgsw_key(key);
  free_trgsw(c_1);
  free_trgsw(c_2);
  free_trgsw(z_inv);
  free_trgsw(one);
}

void test_trgsw_mul_by_xai(){
  TRLWE_Key trlwe_key = trlwe_new_binary_key(N, k, rlwe_std_dev);
  TRGSW_Key key = trgsw_new_key(trlwe_key, l, Bg_bit);
  int test, test2;
  generate_random_bytes(4, (uint8_t *) &test);
  test &= (N - 1);
  generate_random_bytes(4, (uint8_t *) &test2);
  test2 &= (N - 1);
  TRGSW c_1 = trgsw_new_exp_sample(test, key);
  TRGSW res_c = trgsw_alloc_new_sample(l, Bg_bit, k, N);
  trgsw_mul_by_xai(res_c, c_1, test2);
  uint64_t res = _debug_trgsw_decrypt_exp_sample(res_c, key);
  TEST_ASSERT_EQUAL_INT64_MESSAGE((test+test2)&(N - 1), res, "TRGSW MUL by X^a FAILED.");
  free_trlwe_key(trlwe_key);
  free_trgsw_key(key);
  free_trgsw(c_1);
  free_trgsw(res_c);
}

void test_trgsw_dft(){
  TRLWE_Key trlwe_key = trlwe_new_binary_key(N, k, rlwe_std_dev);
  int test;
  generate_random_bytes(4, (uint8_t *) &test);
  test &= (N - 1);
  TRGSW_Key key = trgsw_new_key(trlwe_key, l, Bg_bit);
  TRGSW c_1 = trgsw_new_exp_sample(test, key);
  TRGSW c_2 = trgsw_alloc_new_sample(l, Bg_bit, k, N);
  TRGSW_DFT c_1_dft = trgsw_alloc_new_DFT_sample(l, Bg_bit, k, N);
  trgsw_to_DFT(c_1_dft, c_1);
  trgsw_from_DFT(c_2, c_1_dft);
  trgsw_to_DFT(c_1_dft, c_2);
  trgsw_from_DFT(c_2, c_1_dft);
  trgsw_to_DFT(c_1_dft, c_2);
  trgsw_from_DFT(c_2, c_1_dft);
  trgsw_to_DFT(c_1_dft, c_2);
  trgsw_from_DFT(c_2, c_1_dft);
  uint64_t res = _debug_trgsw_decrypt_exp_sample(c_2, key);
  TEST_ASSERT_EQUAL_INT64_MESSAGE(test, res, "TRGSW DFT FAILED.");
  free_trlwe_key(trlwe_key);
  free_trgsw_key(key);
  free_trgsw(c_1);
  free_trgsw(c_2);
  free_trgsw(c_1_dft);
}

void test_trgsw_trlwe_mul(){
  TRLWE_Key trlwe_key = trlwe_new_binary_key(N, k, rlwe_std_dev);
  int test;
  generate_random_bytes(4, (uint8_t *) &test);
  test &= (N - 1);
  TRGSW_Key key = trgsw_new_key(trlwe_key, l, Bg_bit);
  TRGSW c_1 = trgsw_new_exp_sample(test, key);
  TRGSW_DFT c_1_dft = trgsw_alloc_new_DFT_sample(l, Bg_bit, k, N);
  trgsw_to_DFT(c_1_dft, c_1);
  TorusPolynomial poly_1 = polynomial_new_torus_polynomial(N);
  TorusPolynomial poly_2 = polynomial_new_torus_polynomial(N);
  TorusPolynomial poly_res = polynomial_new_torus_polynomial(N);
  TorusPolynomial poly_clear = polynomial_new_torus_polynomial(N);
  generate_random_bytes(N*sizeof(Torus), (uint8_t *) poly_1->coeffs);
  TRLWE c_2 = trlwe_new_sample(poly_1, trlwe_key);
  TRLWE c_res = trlwe_new_sample(poly_1, trlwe_key);
  TRLWE_DFT c_2_dft = trlwe_alloc_new_DFT_sample(k, N);
  trgsw_mul_trlwe_DFT(c_2_dft, c_2, c_1_dft);
  trlwe_from_DFT(c_res, c_2_dft);
  trlwe_phase(poly_res, c_res, trlwe_key);
  memset(poly_clear->coeffs, 0, N*sizeof(Torus));
  memset(poly_2->coeffs, 0, N*sizeof(Torus));
  poly_2->coeffs[test] = 1;
  polynomial_naive_mul_addto_torus(poly_clear, poly_1, poly_2);
  TEST_ASSERT_TORUS_ARRAY_WITHIN_MESSAGE(1UL << 54, poly_clear->coeffs, poly_res->coeffs, N, "External product failed.");
  free_trlwe_key(trlwe_key);
  free_trgsw_key(key);
  free_trgsw(c_1);
  free_trgsw(c_1_dft);
  free_trlwe(c_2);
  free_trlwe(c_2_dft);
  free_trlwe(c_res);
  free_polynomial(poly_1);
  free_polynomial(poly_2);
  free_polynomial(poly_res);
  free_polynomial(poly_clear);
}

void test_trgsw_mul(){
  TRLWE_Key trlwe_key = trlwe_new_binary_key(N, k, rlwe_std_dev);
  int test[3];
  generate_random_bytes(12, (uint8_t *) test);
  test[0] &= (N - 1);
  test[1] &= (N - 1);
  TRGSW_Key key = trgsw_new_key(trlwe_key, l, Bg_bit);
  TRGSW c_1 = trgsw_new_exp_sample(test[0], key);
  TRGSW c_2 = trgsw_new_exp_sample(test[1], key);
  TRGSW_DFT c_1_dft = trgsw_alloc_new_DFT_sample(l, Bg_bit, k, N);
  TRGSW_DFT c_2_dft = trgsw_alloc_new_DFT_sample(l, Bg_bit, k, N);
  TRGSW_DFT c_3_dft = trgsw_alloc_new_DFT_sample(l, Bg_bit, k, N);
  trgsw_to_DFT(c_1_dft, c_1);
  trgsw_to_DFT(c_2_dft, c_2);
  trgsw_mul_DFT2(c_3_dft, c_1_dft, c_2_dft);
  trgsw_from_DFT(c_1, c_3_dft);
  int res = _debug_trgsw_decrypt_exp_sample(c_1, key);
  TEST_ASSERT_EQUAL_INT64_MESSAGE((test[0]+test[1])&(N - 1), res, "TRGSW MUL FAILED.");
  free_trlwe_key(trlwe_key);
  free_trgsw_key(key);
  free_trgsw(c_1);
  free_trgsw(c_2);
  free_trgsw(c_1_dft);
  free_trgsw(c_2_dft);
  free_trgsw(c_3_dft);
}

void test_trgsw_reg_sub(){
  TRLWE_Key trlwe_key = trlwe_new_binary_key(N, k, rlwe_std_dev);
  TRGSW_Key key = trgsw_new_key(trlwe_key, l, Bg_bit);

  int test[3];
  generate_random_bytes(12, (uint8_t *) test);
  test[0] &= (N - 1);
  test[1] &= (N - 1);

  TRGSW_REG * reg = trgsw_reg_alloc_array(3, l, Bg_bit, k, N);
  TRGSW tmp = trgsw_new_exp_sample(test[0], key);
  trgsw_to_DFT(reg[0]->positive, tmp);
  free_trgsw(tmp);
  tmp = trgsw_new_exp_sample(N - test[0], key);
  trgsw_to_DFT(reg[0]->negative, tmp);
  free_trgsw(tmp);

  tmp = trgsw_new_exp_sample(test[1], key);
  trgsw_to_DFT(reg[1]->positive, tmp);
  free_trgsw(tmp);
  tmp = trgsw_new_exp_sample(N - test[1], key);
  trgsw_to_DFT(reg[1]->negative, tmp);

  
  trgsw_reg_sub(reg[2], reg[0], reg[1]);
  trgsw_from_DFT(tmp, reg[2]->positive);

  int res = _debug_trgsw_decrypt_exp_sample(tmp, key);
  TEST_ASSERT_EQUAL_INT64_MESSAGE((test[0]-test[1])&(N - 1), res, "TRGSW REG SUB FAILED.");

  trgsw_reg_subto(reg[0], reg[1]);
  trgsw_from_DFT(tmp, reg[0]->positive);

  res = _debug_trgsw_decrypt_exp_sample(tmp, key);
  TEST_ASSERT_EQUAL_INT64_MESSAGE((test[0]-test[1])&(N - 1), res, "TRGSW REG SUBTO FAILED.");

  for (size_t i = 0; i < 99; i++){
    trgsw_reg_subto(reg[0], reg[1]);
  }
  
  trgsw_from_DFT(tmp, reg[0]->positive);

  res = _debug_trgsw_decrypt_exp_sample(tmp, key);
  TEST_ASSERT_EQUAL_INT64_MESSAGE((test[0]-100*test[1])&(N - 1), res, "TRGSW REG 1000x SUBTO FAILED.");

  int result = test[1];
  for (size_t i = 0; i < 10; i++){
    generate_random_bytes(4, (uint8_t *) test);
    test[0] &= (N - 1);
    trgsw_monomial_sample(tmp, 1, test[0], key);
    trgsw_to_DFT(reg[0]->positive, tmp);
    trgsw_monomial_sample(tmp, 1, N - test[0], key);
    trgsw_to_DFT(reg[0]->negative, tmp);
    
    trgsw_reg_subto(reg[1], reg[0]);
    result = (result - test[0]) & (N - 1);
  }
  
  trgsw_from_DFT(tmp, reg[1]->positive);

  res = _debug_trgsw_decrypt_exp_sample(tmp, key);
  TEST_ASSERT_EQUAL_INT64_MESSAGE(result, res, "TRGSW REG 10x Random SUBTO FAILED.");
  
  free_trlwe_key(trlwe_key);
  free_trgsw_key(key);
  free_trgsw(tmp);
  free_trgsw_reg_array(reg, 3);
}

void test_trgsw_reg_sub2(){
  TRLWE_Key trlwe_key = trlwe_new_binary_key(N, k, rlwe_std_dev);
  TRGSW_Key key = trgsw_new_key(trlwe_key, l, Bg_bit);

  uint64_t a[6], z_c[6], test[11];
  const int size = 6;
  generate_random_bytes(6*8, (uint8_t *) a);
  generate_random_bytes(6*8, (uint8_t *) z_c);

  TRGSW z[6];

  for (size_t i = 0; i < size; i++){
    a[i] &= (N-1);
    z_c[i] &= 1;
    z[i] = trgsw_new_sample(z_c[i], key);
  }
  

  TRGSW tmp = trgsw_alloc_new_sample(l, Bg_bit, k, N), z_inv = trgsw_alloc_new_sample(l, Bg_bit, k, N), one = trgsw_new_noiseless_trivial_sample(1, l, Bg_bit, k, N);
  TRGSW_REG * out = trgsw_reg_alloc_array(11, l, Bg_bit, k, N);
  TRGSW_REG reg_tmp = trgsw_reg_alloc(l, Bg_bit, k, N);


  for (size_t i = 0; i < size; i++){
    for (size_t j = 0; j < size; j++){
      trgsw_sub(z_inv, one, z[i]);
      trgsw_mul_by_xai(tmp, z[i], a[j]);
      trgsw_addto(tmp, z_inv);
      if(i == 0 || j == size - 1){
        test[i + j] = z_c[i]*a[j];
        trgsw_to_DFT(out[i+j]->positive, tmp);
      }else{
        trgsw_to_DFT(reg_tmp->positive, tmp);
      }
      trgsw_mul_by_xai(tmp, z[i], N - a[j]);
      trgsw_addto(tmp, z_inv);
      if(i == 0 || j == size - 1){
        trgsw_to_DFT(out[i+j]->negative, tmp);
      }else{
        test[i + j] += z_c[i]*a[j];
        trgsw_to_DFT(reg_tmp->negative, tmp);
        trgsw_reg_addto(out[i+j], reg_tmp);
      }
      trgsw_from_DFT(tmp, out[i+j]->positive);
    }
  }



  int d = 9;
  // for (int i = 10; i >= 6; i--){
  //   trgsw_reg_subto(out[i - d/3], out[i]);
  //   trgsw_reg_subto(out[i - 2*d/3], out[i]);
  //   test[i-d/3] = (test[i-d/3] - test[i])&(N - 1);
  //   test[i-2*d/3] = (test[i-2*d/3] - test[i])&(N - 1);
  //   printf("[%d] %lu %lu \n", i, test[i-2*d/3], _debug_trgsw_decrypt_exp_DFT_sample(out[i - 2*d/3]->positive, key));
  // }



  trgsw_reg_subto(out[9 - d/3], out[9]);
  test[9-d/3] = (test[9-d/3] - test[9])&(N - 1);

  trgsw_reg_sample(reg_tmp, 56, key);
  trgsw_reg_subto(reg_tmp, out[6]);

  TEST_ASSERT_EQUAL_HEX64_MESSAGE((56 - test[6])&(N - 1), _debug_trgsw_decrypt_exp_DFT_sample(reg_tmp->positive, key), "TRGSW REG Single Sequential SUBTO FAILED.");
  return;
  uint64_t res[11];
  for (int i = 0; i < 11; i++){
    res[i] = _debug_trgsw_decrypt_exp_DFT_sample(out[i]->positive, key);
  }
  TEST_ASSERT_EQUAL_HEX64_ARRAY_MESSAGE(test, res, 11, "TRGSW REG Sequential SUBTO FAILED.");
  free_trlwe_key(trlwe_key);
  free_trgsw_key(key);
  free_trgsw(tmp);
  free_trgsw_reg(reg_tmp);
  free_trgsw_reg_array(out, 11);
  for (size_t i = 0; i < size; i++) free_trgsw(z[i]);
}

void print_TRGSW_DFT(TRGSW_DFT p, int size, TRGSW_Key key){
  int first_print = 1;
  const uint64_t N = key->trlwe_key->s[0]->N, k = key->trlwe_key->k;
  TRLWE tmp = trlwe_new_noiseless_trivial_sample(NULL, k, N);
  tmp->b->coeffs[0] = (1UL << (sizeof(Torus)*8 - key->Bg_bit));
  TRLWE_DFT res = trlwe_alloc_new_DFT_sample(k, N);
  trgsw_mul_trlwe_DFT(res, tmp, p);
  trlwe_from_DFT(tmp, res);
  TorusPolynomial poly = polynomial_new_torus_polynomial(N);
  trlwe_phase(poly, tmp, key->trlwe_key);


  for (int i = size - 1; i >=0; i--){
    if(!first_print) printf(" + ");
    else first_print = 0;
    uint64_t res = poly->coeffs[i];
    // if(res == -1){
    //   printf("TRGSW REG Polynomial dec Failed. No value.\n");
    //   return;
    // }
    printf("%lx*x^%d", res, i);
  }
  printf("\n");
  free_polynomial(poly);
}

void test_trgsw_reg_sub3(){
  TRLWE_Key trlwe_key = trlwe_new_binary_key(N, k, rlwe_std_dev);
  TRGSW_Key key = trgsw_new_key(trlwe_key, l, Bg_bit);

  const int size = 10;
  TRGSW_DFT * c = trgsw_alloc_new_DFT_sample_array(size, l, Bg_bit, k, N);
  TRGSW tmp = trgsw_alloc_new_sample(l, Bg_bit, k, N);

  uint64_t test[size];
  generate_random_bytes(size*8, (uint8_t *)test);
  
  for (size_t i = 0; i < size; i++){
    test[i] &= (N - 1);
    trgsw_monomial_sample(tmp, 1, test[i], key);
    trgsw_to_DFT(c[i], tmp);
  }

  trgsw_mul_DFT2(c[4], c[1], c[0]);
  trgsw_mul_DFT2(c[5], c[2], c[4]);
  trgsw_mul_DFT2(c[6], c[3], c[5]);


  // print_TRGSW_DFT(c[4], N, key);
  // print_TRGSW_DFT(c[5], N, key);
  // print_TRGSW_DFT(c[6], N, key);

  
  TEST_ASSERT_EQUAL_HEX64_MESSAGE((test[0] + test[1])&(N - 1), _debug_trgsw_decrypt_exp_DFT_sample(c[4], key), "TRGSW Sequential SUBTO FAILED 1");
  TEST_ASSERT_EQUAL_HEX64_MESSAGE((test[0] + test[1] + test[2])&(N - 1), _debug_trgsw_decrypt_exp_DFT_sample(c[5], key), "TRGSW Sequential SUBTO FAILED 2");
  TEST_ASSERT_EQUAL_HEX64_MESSAGE((test[0] + test[1] + test[2] + test[3])&(N - 1), _debug_trgsw_decrypt_exp_DFT_sample(c[6], key), "TRGSW Sequential SUBTO FAILED 3");
  free_trlwe_key(trlwe_key);
  free_trgsw_key(key);
  free_trgsw(tmp);
  free_trgsw_array(c, size);
}


void print_TRGSW(TRGSW p, int size, TRGSW_Key key){
  int first_print = 1;
  const uint64_t N = key->trlwe_key->s[0]->N, k = key->trlwe_key->k;
  TRLWE tmp = trlwe_new_noiseless_trivial_sample(NULL, k, N);
  tmp->b->coeffs[0] = (1UL << (sizeof(Torus)*8 - key->Bg_bit));
  TRLWE res = trlwe_alloc_new_sample(k, N);
  trgsw_naive_mul_trlwe(res, tmp, p);
  TorusPolynomial poly = polynomial_new_torus_polynomial(N);
  trlwe_phase(poly, res, key->trlwe_key);


  for (int i = size - 1; i >=0; i--){
    if(!first_print) printf(" + ");
    else first_print = 0;
    uint64_t res = poly->coeffs[i];
    // if(res == -1){
    //   printf("TRGSW REG Polynomial dec Failed. No value.\n");
    //   return;
    // }
    printf("%lx*x^%d", res, i);
  }
  printf("\n");
  free_polynomial(poly);
}


void test_trgsw_naive_mul(){
  TRLWE_Key trlwe_key = trlwe_new_binary_key(N, k, rlwe_std_dev);
  TRGSW_Key key = trgsw_new_key(trlwe_key, l, Bg_bit);

  const int size = 10;
  TRGSW * c = trgsw_alloc_new_sample_array(size, l, Bg_bit, k, N);

  uint64_t test[size];
  generate_random_bytes(size*8, (uint8_t *)test);
  
  for (size_t i = 0; i < size; i++){
    test[i] &= (N - 1);
    trgsw_monomial_sample(c[i], 1, test[i], key);
  }

  trgsw_naive_mul(c[4], c[1], c[0]);
  // print_TRGSW(c[4], N, key);

  trgsw_naive_mul(c[5], c[2], c[4]);
  // print_TRGSW(c[5], N, key);

  trgsw_naive_mul(c[6], c[3], c[5]);
  // print_TRGSW(c[6], N, key);

  
  TEST_ASSERT_EQUAL_HEX64_MESSAGE((test[0] + test[1])&(N - 1), _debug_trgsw_decrypt_exp_sample(c[4], key), "TRGSW Sequential SUBTO FAILED 1");
  TEST_ASSERT_EQUAL_HEX64_MESSAGE((test[0] + test[1] + test[2])&(N - 1), _debug_trgsw_decrypt_exp_sample(c[5], key), "TRGSW Sequential SUBTO FAILED 2");
  TEST_ASSERT_EQUAL_HEX64_MESSAGE((test[0] + test[1] + test[2] + test[3])&(N - 1), _debug_trgsw_decrypt_exp_sample(c[6], key), "TRGSW Sequential SUBTO FAILED 3");
  free_trlwe_key(trlwe_key);
  free_trgsw_key(key);
  free_trgsw_array(c, size);
}


void _test_tlwe_ks_save_keys_and_exit(){
  TLWE_Key key_1 = tlwe_new_binary_key(n, lwe_std_dev);
  TLWE_Key key_2 = tlwe_new_binary_key(N, rlwe_std_dev);
  TLWE_KS_Key ks_key_2_to_1 = tlwe_new_KS_key(key_1, key_2, t, base_bit);
  FILE * f_test = fopen("test_key", "w+");
  tlwe_save_key(f_test, key_1);
  tlwe_save_key(f_test, key_2);
  tlwe_save_KS_key(f_test, ks_key_2_to_1);
  fclose(f_test);
  exit(0);
}


void test_tlwe_ks(){
  TLWE_Key key_1 = tlwe_new_binary_key(n, lwe_std_dev);
  TLWE_Key key_2 = tlwe_new_binary_key(N, rlwe_std_dev);
  TLWE_KS_Key ks_key_2_to_1 = tlwe_new_KS_key(key_1, key_2, t, base_bit);
  // _test_tlwe_ks_save_keys_and_exit();

  // FILE * f_test = fopen("test_key", "r");
  // TLWE_Key key_1 = tlwe_load_new_key(f_test);
  // TLWE_Key key_2 = tlwe_load_new_key(f_test);
  // // TLWE_KS_Key ks_key_2_to_1 = tlwe_load_new_KS_key(f_test);
  // fclose(f_test);
  // TLWE_KS_Key ks_key_2_to_1 = tlwe_new_KS_key(key_1, key_2, t, base_bit);


  Torus in1, in2;
  generate_random_bytes(sizeof(Torus), (uint8_t *) &in2);
  TLWE c_2 = tlwe_new_sample(in2, key_2);
  TLWE c_2_k1 = tlwe_new_noiseless_trivial_sample(0, key_1->n);
  // printf("in2: %lx, c2: %lx c_2_k1: %lx\n", in2, tlwe_phase(c_2, key_2), tlwe_phase(c_2_k1, key_1));
  tlwe_keyswitch(c_2_k1, c_2, ks_key_2_to_1);
  // printf("in2: %lx, c2: %lx c_2_k1: %lx\n", in2, tlwe_phase(c_2, key_2), tlwe_phase(c_2_k1, key_1));
  TEST_ASSERT_TORUS_WITHIN_MESSAGE(1UL << 58, in2, tlwe_phase(c_2_k1, key_1), "TLWE ks_key_2_to_1 failed");

  TLWE_KS_Key ks_key_1_to_2 = tlwe_new_KS_key(key_2, key_1, t, base_bit);
  generate_random_bytes(sizeof(Torus), (uint8_t *) &in1);
  TLWE c_1 = tlwe_new_sample(in1, key_1);
  TLWE c_1_k2 = tlwe_new_noiseless_trivial_sample(0, key_2->n);
  // printf("in1: %lx, c1: %lx\n", in1, tlwe_phase(c_1, key_1));
  tlwe_keyswitch(c_1_k2, c_1, ks_key_1_to_2);
  TEST_ASSERT_TORUS_WITHIN_MESSAGE(1UL << 58, in1, tlwe_phase(c_1_k2, key_2), "TLWE ks_key_1_to_2 failed");

  free_tlwe(c_1);
  free_tlwe(c_2);
  free_tlwe(c_1_k2);
  free_tlwe(c_2_k1);
  free_tlwe_ks_key(ks_key_1_to_2);
  free_tlwe_ks_key(ks_key_2_to_1);
  free_tlwe_key(key_1);
  free_tlwe_key(key_2);
}



void test_tlwe_pack1_ks(){
  TLWE_Key key_1 = tlwe_new_binary_key(n, lwe_std_dev);
  TRLWE_Key key_2 = trlwe_new_binary_key(N, k, rlwe_std_dev);
  Generic_KS_Key ksk = trlwe_new_packing1_KS_key(key_2, key_1, t, base_bit);

  Torus in1;
  generate_random_bytes(sizeof(Torus), (uint8_t *) &in1);
  TLWE c_1 = tlwe_new_sample(in1, key_1);
  TRLWE out_enc = trlwe_alloc_new_sample(k, N); 

  trlwe_packing1_keyswitch(out_enc, c_1, ksk);

  TorusPolynomial out = polynomial_new_torus_polynomial(N);
  TorusPolynomial res = polynomial_new_torus_polynomial(N);
  memset(res->coeffs, 0, N*sizeof(Torus));
  res->coeffs[0] = in1;

  trlwe_phase(out, out_enc, key_2);

  TEST_ASSERT_TORUS_ARRAY_WITHIN_MESSAGE(1UL << 58, res->coeffs, out->coeffs, N, "TLWE(m) -> TRLWE(m*X^0) Key Switching failed");

  free_polynomial(out);
  free_polynomial(res);
  free_tlwe(c_1);
  free_tlwe_key(key_1);
  free_trlwe(out_enc);
  free_trlwe_key(key_2);
  free_trlwe_generic_ks_key(ksk);
}


void test_tlwe_pack1_ks_CDKS21(){
  TLWE_Key key_1 = tlwe_new_binary_key(N, rlwe_std_dev);
  TRLWE_Key key_2 = trlwe_new_binary_key(N, k, rlwe_std_dev);
  trlwe_extract_tlwe_key(key_1, key_2);
  TRLWE_KS_Key * ksk = trlwe_new_packing1_KS_key_CDKS21(key_2, key_1, l, Bg_bit);
  _glb_debug_trlwe_key = key_2;

  Torus in1;
  generate_random_bytes(sizeof(Torus), (uint8_t *) &in1);
  // printf("in1: %lx\n", in1);

  TLWE c_1 = tlwe_new_sample(in1, key_1);
  TRLWE out_enc = trlwe_alloc_new_sample(k, N); 

  trlwe_packing1_keyswitch_CDKS21(out_enc, c_1, ksk);

  TorusPolynomial out = polynomial_new_torus_polynomial(N);
  TorusPolynomial res = polynomial_new_torus_polynomial(N);
  memset(res->coeffs, 0, N*sizeof(Torus));
  res->coeffs[0] = in1*N;

  trlwe_phase(out, out_enc, key_2);

  // _debug_print_TRLWE(out_enc, key_2);

  TEST_ASSERT_TORUS_ARRAY_WITHIN_MESSAGE(1UL << 58, res->coeffs, out->coeffs, N, "TLWE(m) -> TRLWE(m*X^0) Key Switching failed");

  free_polynomial(out);
  free_polynomial(res);
  free_tlwe(c_1);
  free_tlwe_key(key_1);
  free_trlwe(out_enc);
  free_trlwe_key(key_2);
  for (size_t i = 0; i < log2(N); i++){
    free_trlwe_ks_key(ksk[i]);
  }
}

void test_tlwe_pack_key_priv_ks(){
  const int bit_len = sizeof(Torus)*8;
  TLWE_Key key_1 = tlwe_new_binary_key(n, lwe_std_dev);
  TRLWE_Key key_2 = trlwe_new_binary_key(N, k, rlwe_std_dev);
  Generic_KS_Key ksk = trlwe_new_priv_SK_KS_key(key_2, key_1, t, base_bit);

  Torus in1 = 1UL << (bit_len - Bg_bit);
  // generate_random_bytes(sizeof(Torus), (uint8_t *) &in1);
  TLWE c_1 = tlwe_new_sample(in1, key_1);
  TRLWE out_enc = trlwe_alloc_new_sample(k, N); 

  trlwe_priv_keyswitch(out_enc, c_1, ksk);

  TorusPolynomial out = polynomial_new_torus_polynomial(N);
  TorusPolynomial res = polynomial_new_torus_polynomial(N);
  memset(res->coeffs, 0, N);
  for (size_t i = 0; i < N; i++){
    res->coeffs[i] = in1 * -key_2->s[0]->coeffs[i];
  }
  
  trlwe_phase(out, out_enc, key_2);

  TEST_ASSERT_TORUS_ARRAY_WITHIN_MESSAGE(1UL << 58, res->coeffs, out->coeffs, N, "Private Key Switching TLWE(M) -> TRLWE(m*-s) failed");

  free_polynomial(out);
  free_polynomial(res);
  free_tlwe(c_1);
  free_tlwe_key(key_1);
  free_trlwe(out_enc);
  free_trlwe_key(key_2);
  free_trlwe_generic_ks_key(ksk);
}

void test_trlwe_pack_key_priv_ks(){
  const int bit_len = sizeof(Torus)*8;
  TRLWE_Key key = trlwe_new_binary_key(N, k, rlwe_std_dev);
  TRLWE_KS_Key * ksk = trlwe_new_priv_KS_key(key, key, t, base_bit);

  Torus in1 = 1UL << (bit_len - Bg_bit);
  // generate_random_bytes(sizeof(Torus), (uint8_t *) &in1);
  TRLWE c_1 = trlwe_new_sample(NULL, key);
  c_1->b->coeffs[0] += in1;
  
  TRLWE out_enc = trlwe_alloc_new_sample(k, N); 

  trlwe_priv_keyswitch_2(out_enc, c_1, ksk);

  TorusPolynomial out = polynomial_new_torus_polynomial(N);
  TorusPolynomial res = polynomial_new_torus_polynomial(N);
  memset(res->coeffs, 0, N);
  for (size_t i = 0; i < N; i++){
    res->coeffs[i] = in1 * -key->s[0]->coeffs[i];
  }
  
  trlwe_phase(out, out_enc, key);

  TEST_ASSERT_TORUS_ARRAY_WITHIN_MESSAGE(1UL << 52, res->coeffs, out->coeffs, N, "Private Key Switching TRLWE(M) -> TRLWE(m*-s) failed");

  free_polynomial(out);
  free_polynomial(res);
  free_trlwe(c_1);
  free_trlwe(out_enc);
  free_trlwe_key(key);
  free_trlwe_ks_key(ksk[0]);
  free_trlwe_ks_key(ksk[1]);

}

void test_multivalue_bootstrap_CLOT21(){
  TLWE_Key key_tlwe = tlwe_new_binary_key(n, lwe_std_dev);
  TLWE_Key key_tlwe_out = tlwe_new_binary_key(N, rlwe_std_dev);
  TRLWE_Key key_trlwe = trlwe_new_binary_key(N, k, rlwe_std_dev);
  trlwe_extract_tlwe_key(key_tlwe_out, key_trlwe);
  _glb_debug_trlwe_key = key_trlwe;
  _glb_debug_tlwe_key = key_tlwe_out;

  TRGSW_Key trgsw_key = trgsw_new_key(key_trlwe, l, Bg_bit);
  Bootstrap_Key bk_key = new_bootstrap_key(trgsw_key, key_tlwe, 1);

  TLWE * out = tlwe_alloc_sample_array(8, N);
  TRLWE tv = trlwe_alloc_new_sample(k, N);
  Torus lut[8*2];
  generate_random_bytes(sizeof(Torus)*2*8, (uint8_t *) lut);
  trlwe_torus_packing(tv, lut, 2*8);
  
  TLWE sel = tlwe_new_sample(double2torus(1./4), key_tlwe);
  multivalue_bootstrap_CLOT21(out, tv, sel, bk_key, 2, 8);

  for (size_t i = 0; i < 8; i++){
    TEST_ASSERT_TORUS_WITHIN_MESSAGE(1UL << 58, lut[8 + i], tlwe_phase(out[i], key_tlwe_out), "Multi value bootstrap CLOT21 failed.");
  }
  fflush(stdout);
  free_tlwe_key(key_tlwe);
  free_tlwe_key(key_tlwe_out);
  free_trlwe_key(key_trlwe);
  free_trgsw_key(trgsw_key);
  free_bootstrap_key(bk_key);
  free_trlwe(tv);
  free_tlwe(sel);
  free_tlwe_array(out, 8);
}

void test_circuit_bootstrap(){
  SKIP_IF_TORUS32
  const int l = 4, Bg_bit = 9, t = 6, base_bit = 4; // Warning: High memory consuming 
  TLWE_Key key_1 = tlwe_new_binary_key(n, lwe_std_dev);
  TRLWE_Key key_2 = trlwe_new_binary_key(N, k, rlwe_std_dev);
  TRGSW_Key key_3 = trgsw_new_key(key_2, l, Bg_bit);
  TLWE_Key key_2_extracted = tlwe_new_binary_key(N, rlwe_std_dev);
  trlwe_extract_tlwe_key(key_2_extracted, key_2);
  TRLWE_KS_Key * kska = trlwe_new_priv_KS_key(key_2, key_2, 20, 2);
  Generic_KS_Key kskb = trlwe_new_packing1_KS_key(key_2, key_2_extracted, t, base_bit);
  TRGSW out = trgsw_alloc_new_sample(l, Bg_bit, k, N);
  Bootstrap_Key bk = new_bootstrap_key(key_3, key_1, 1);
  TorusPolynomial * p = polynomial_new_array_of_torus_polynomials(N, 2); 

  // Test LWE(1/4) -> TRGSW (1)
  TLWE c_in = tlwe_new_sample(double2torus(1./4), key_1);
  circuit_bootstrap_3(out, c_in, bk, kska, kskb);


  generate_random_bytes(N*sizeof(Torus), (uint8_t *) p[0]->coeffs);
  // TRLWE rnd = trlwe_new_noiseless_trivial_sample(p[0], k, N);
  TRLWE rnd = trlwe_new_sample(p[0], key_2);
  TRGSW_DFT out_dft = trgsw_alloc_new_DFT_sample(l, Bg_bit, k, N);
  TRLWE_DFT mul_res_dft = trlwe_alloc_new_DFT_sample(k, N);
  TRLWE mul_res = trlwe_alloc_new_sample(k, N);
  trgsw_to_DFT(out_dft, out);
  trgsw_mul_trlwe_DFT(mul_res_dft, rnd, out_dft);
  trlwe_from_DFT(mul_res, mul_res_dft);
  trlwe_phase(p[1], mul_res, key_2);
  TEST_ASSERT_TORUS_ARRAY_WITHIN_MESSAGE(1UL << 58, p[0]->coeffs, p[1]->coeffs, N, "Circuit Bootstrap LWE(1/4) -> TRGSW (1) failed");
  
  // Test LWE(0) -> TRGSW (0)
  tlwe_sample(c_in, 0, key_1);

  circuit_bootstrap_3(out, c_in, bk, kska, kskb);

  generate_random_bytes(N*sizeof(Torus), (uint8_t *) p[0]->coeffs);
  trgsw_to_DFT(out_dft, out);
  trgsw_mul_trlwe_DFT(mul_res_dft, rnd, out_dft);
  trlwe_from_DFT(mul_res, mul_res_dft);
  trlwe_phase(p[1], mul_res, key_2);
  memset(p[0]->coeffs, 0, N*sizeof(Torus));
  TEST_ASSERT_TORUS_ARRAY_WITHIN_MESSAGE(1UL << 58, p[0]->coeffs, p[1]->coeffs, N, "Circuit Bootstrap LWE(0/4) -> TRGSW (0) failed");
  
  free_array_of_polynomials((void **) p, 2);
  free_trgsw(out_dft);
  free_trlwe(rnd);
  free_trlwe(mul_res);
  free_trlwe(mul_res_dft);
  free_trgsw(out);
  free_tlwe_key(key_1);
  free_trlwe_key(key_2);
  free_trgsw_key(key_3);
  free_trlwe_ks_key(kska[0]);
  free_trlwe_ks_key(kska[1]);
  free_trlwe_generic_ks_key(kskb);
  free_tlwe(c_in);
}

void test_public_mux(){
  const int bit_len = sizeof(Torus)*8;
  TRLWE_Key key = trlwe_new_binary_key(N, k, 0);
  TorusPolynomial * p = polynomial_new_array_of_torus_polynomials(N, 3);
  TRLWE out = trlwe_alloc_new_sample(k, N);
  TRLWE_DFT * selector = trlwe_alloc_new_DFT_sample_array(l, k, N);
  TRLWE tmp_trlwe = trlwe_alloc_new_sample(k, N);
  for (size_t i = 0; i < l; i++){
    Torus sign = 1UL << (bit_len - (i + 1) * Bg_bit);
    trlwe_sample(tmp_trlwe, NULL, key);
    tmp_trlwe->b->coeffs[0] += sign;
    trlwe_to_DFT(selector[i], tmp_trlwe);
  }
  public_mux(out, p[0], p[1], selector, l, Bg_bit);
  trlwe_phase(p[2], out, key);
  TEST_ASSERT_TORUS_ARRAY_WITHIN_MESSAGE(1UL << 56, p[1]->coeffs, p[2]->coeffs, N, "Public MUX s=1 failed");

  for (size_t i = 0; i < l; i++){
    Torus sign = 0UL << (bit_len - (i + 1) * Bg_bit);
    trlwe_sample(tmp_trlwe, NULL, key);
    tmp_trlwe->b->coeffs[0] += sign;
    trlwe_to_DFT(selector[i], tmp_trlwe);
  }
  public_mux(out, p[0], p[1], selector, l, Bg_bit);
  trlwe_phase(p[2], out, key);
  TEST_ASSERT_TORUS_ARRAY_WITHIN_MESSAGE(1UL << 56, p[0]->coeffs, p[2]->coeffs, N, "Public MUX s=0 failed");

  free_trlwe_key(key);
  free_array_of_polynomials((void **) p, 3);
  free_trlwe(out);
  free_trlwe(tmp_trlwe);
  free_trlwe_array(selector, l);
}

void test_FDFB_KS21(){
  const int t = 6, base_bit = 4;
  TLWE_Key key_tlwe = tlwe_new_binary_key(n, lwe_std_dev);
  TLWE_Key key_tlwe_out = tlwe_new_binary_key(N, rlwe_std_dev);
  TRLWE_Key key_trlwe = trlwe_new_binary_key(N, k, rlwe_std_dev);
  trlwe_extract_tlwe_key(key_tlwe_out, key_trlwe);
  Generic_KS_Key ksk = trlwe_new_packing1_KS_key(key_trlwe, key_tlwe_out, t, base_bit);
  _glb_debug_trlwe_key = key_trlwe;

  TRGSW_Key trgsw_key = trgsw_new_key(key_trlwe, l, Bg_bit);
  Bootstrap_Key bk_key = new_bootstrap_key(trgsw_key, key_tlwe, 1);

  Torus in[8];
  TLWE c_in = tlwe_alloc_sample(n);
  TLWE c_out = tlwe_alloc_sample(N);
  generate_random_bytes(sizeof(Torus)*8, (uint8_t *) in);

  TorusPolynomial poly_res = polynomial_new_torus_polynomial(2*N);
  for (size_t i = 0; i < 2*N; i++) poly_res->coeffs[i] = in[i/(N/4)];

  for (size_t i = 0; i < 8; i++){
    tlwe_sample(c_in, int2torus(i, 3), key_tlwe);
    full_domain_functional_bootstrap_KS21(c_out, poly_res, c_in, bk_key, ksk, 8);
    TEST_ASSERT_TORUS_WITHIN_MESSAGE(1UL << 58, in[i], tlwe_phase(c_out, key_tlwe_out), "FDFB [KS21] failed.");
  }
  
  free_tlwe_key(key_tlwe);
  free_tlwe_key(key_tlwe_out);
  free_trlwe_key(key_trlwe);
  free_trlwe_generic_ks_key(ksk);
  free_trgsw_key(trgsw_key);
  free_bootstrap_key(bk_key);
  free_tlwe(c_in);
  free_tlwe(c_out);
  free_polynomial(poly_res);
}

void test_FDFB_new(){
  TLWE_Key key_tlwe = tlwe_new_binary_key(n, lwe_std_dev);
  TLWE_Key key_tlwe_out = tlwe_new_binary_key(N, rlwe_std_dev);
  TRLWE_Key key_trlwe = trlwe_new_binary_key(N, k, rlwe_std_dev);
  trlwe_extract_tlwe_key(key_tlwe_out, key_trlwe);
  TLWE_KS_Key tlwe_ksk = tlwe_new_KS_key(key_tlwe, key_tlwe_out, t, base_bit); 

  TRGSW_Key trgsw_key = trgsw_new_key(key_trlwe, l, Bg_bit);
  Bootstrap_Key bk_key = new_bootstrap_key(trgsw_key, key_tlwe, 1);

  Torus in[8];
  TLWE c_in = tlwe_alloc_sample(n);
  TLWE c_out = tlwe_alloc_sample(N);
  generate_random_bytes(sizeof(Torus)*8, (uint8_t *) in);

  TRLWE tv = trlwe_alloc_new_sample(k, N);
  trlwe_torus_packing_many_LUT(tv, in, 4, 2);

  for (size_t i = 0; i < 8; i++){
    tlwe_sample(c_in, int2torus(i, 3), key_tlwe);
    full_domain_functional_bootstrap(c_out, tv, c_in, bk_key, tlwe_ksk, 3);
    TEST_ASSERT_TORUS_WITHIN_MESSAGE(1UL << 56, in[i], tlwe_phase(c_out, key_tlwe_out), "FDFB failed.");
  }
  
  free_tlwe_key(key_tlwe);
  free_tlwe_key(key_tlwe_out);
  free_trlwe_key(key_trlwe);
  free_trgsw_key(trgsw_key);
  free_bootstrap_key(bk_key);
  free_tlwe(c_in);
  free_tlwe(c_out);
  free_trlwe(tv);
}

void test_FDFB_CLOT21(){
  // const int l = 8, Bg_bit = 8;
  // const double rlwe_std_dev = pow(2, -54);
  TLWE_Key key_tlwe = tlwe_new_binary_key(n, lwe_std_dev);
  TRLWE_Key key_trlwe = trlwe_new_binary_key(N, k, rlwe_std_dev);
  TRLWE_KS_Key rlk = trlwe_new_RL_key(key_trlwe, 2, 20);
  TLWE_Key key_tlwe_out = tlwe_new_binary_key(N, rlwe_std_dev);
  trlwe_extract_tlwe_key(key_tlwe_out, key_trlwe);
  _glb_debug_trlwe_key = key_trlwe;
  _glb_debug_tlwe_key = key_tlwe_out;
  Generic_KS_Key ksk = trlwe_new_packing1_KS_key(key_trlwe, key_tlwe_out, t, base_bit);

  TRGSW_Key trgsw_key = trgsw_new_key(key_trlwe, l, Bg_bit);
  Bootstrap_Key bk_key = new_bootstrap_key(trgsw_key, key_tlwe, 1);

  Torus in[8];
  TLWE c_in = tlwe_alloc_sample(n);
  TLWE c_out = tlwe_alloc_sample(N);
  generate_random_bytes(sizeof(Torus)*8, (uint8_t *) in);
  const int precision = 4, mod_mask = (1<<precision) - 1;

  for (size_t i = 0; i < 8; i++){
    in[i] = int2torus(in[i]&mod_mask, precision);
  }
  
  
  TRLWE tv[2] = {trlwe_new_noiseless_trivial_sample(NULL, k, N),
                 trlwe_new_noiseless_trivial_sample(NULL, k, N)};
  for (size_t i = 0; i < N; i++) tv[0]->b->coeffs[i] = in[i/(N/4)];
  for (size_t i = N; i < 2*N; i++) tv[1]->b->coeffs[i - N] = in[i/(N/4)];

  for (size_t i = 0; i < 8; i++){
    tlwe_sample(c_in, int2torus(i, 3), key_tlwe);
    full_domain_functional_bootstrap_CLOT21(c_out, tv, c_in, bk_key, ksk, rlk, precision);
    // printf("i: %lu, expected: %lx, was: %lx\n", i, in[i], tlwe_phase(c_out, key_tlwe_out));
    TEST_ASSERT_TORUS_WITHIN_MESSAGE(1ULL << (64 - precision - 1), in[i], tlwe_phase(c_out, key_tlwe_out), "FDFB [CLOT21] failed.");
  }
  
  free_tlwe_key(key_tlwe);
  free_tlwe_key(key_tlwe_out);
  free_trlwe_key(key_trlwe);
  free_trlwe_generic_ks_key(ksk);
  free_trgsw_key(trgsw_key);
  free_bootstrap_key(bk_key);
  free_tlwe(c_in);
  free_tlwe(c_out);
  free_trlwe(tv[0]);
  free_trlwe(tv[1]);
}

void test_FDFB_CLOT21_2(){
  // const int l = 8, Bg_bit = 8;
  // const double rlwe_std_dev = pow(2, -54);
  TLWE_Key key_tlwe = tlwe_new_binary_key(n, lwe_std_dev);
  TRLWE_Key key_trlwe = trlwe_new_binary_key(N, k, rlwe_std_dev);
  TRLWE_KS_Key rlk = trlwe_new_RL_key(key_trlwe, 2, 20);
  TLWE_Key key_tlwe_out = tlwe_new_binary_key(N, rlwe_std_dev);
  trlwe_extract_tlwe_key(key_tlwe_out, key_trlwe);
  _glb_debug_trlwe_key = key_trlwe;
  _glb_debug_tlwe_key = key_tlwe_out;
  Generic_KS_Key ksk = trlwe_new_packing1_KS_key(key_trlwe, key_tlwe_out, t, base_bit);

  TRGSW_Key trgsw_key = trgsw_new_key(key_trlwe, l, Bg_bit);
  Bootstrap_Key bk_key = new_bootstrap_key(trgsw_key, key_tlwe, 1);

  Torus in[8];
  TLWE c_in = tlwe_alloc_sample(n);
  TLWE c_out = tlwe_alloc_sample(N);
  generate_random_bytes(sizeof(Torus)*8, (uint8_t *) in);
  const int precision = 4, mod_mask = (1<<precision) - 1;

  for (size_t i = 0; i < 8; i++){
    in[i] = int2torus(in[i]&mod_mask, precision);
  }

  for (size_t i = 0; i < 8; i++){
    tlwe_sample(c_in, int2torus(i, 3), key_tlwe);
    full_domain_functional_bootstrap_CLOT21_2(c_out, in, c_in, bk_key, ksk, rlk, precision);
    // printf("i: %lu, expected: %lx, was: %lx\n", i, in[i], tlwe_phase(c_out, key_tlwe_out));
    TEST_ASSERT_TORUS_WITHIN_MESSAGE(1ULL << (64 - precision - 1), in[i], tlwe_phase(c_out, key_tlwe_out), "FDFB [CLOT21] failed.");
  }
  
  free_tlwe_key(key_tlwe);
  free_tlwe_key(key_tlwe_out);
  free_trlwe_key(key_trlwe);
  free_trlwe_generic_ks_key(ksk);
  free_trgsw_key(trgsw_key);
  free_bootstrap_key(bk_key);
  free_tlwe(c_in);
  free_tlwe(c_out);
}


void test_trlwe_poly_mul(){
  TRLWE_Key key = trlwe_new_binary_key(N, k, rlwe_std_dev);
  TorusPolynomial in1 = polynomial_new_torus_polynomial(N);
  TorusPolynomial in2 = polynomial_new_torus_polynomial(N);
  TorusPolynomial res = polynomial_new_torus_polynomial(N);
  TorusPolynomial res_e = polynomial_new_torus_polynomial(N);
  generate_random_bytes(sizeof(Torus)*N, (uint8_t *) in1->coeffs);
  generate_random_bytes(sizeof(Torus)*N, (uint8_t *) in2->coeffs);
  for (size_t i = 0; i < N; i++){
    in2->coeffs[i] &= 0xF;
  }
  DFT_Polynomial in2_dft = polynomial_new_DFT_polynomial(N);
  polynomial_naive_mul_torus(res_e, in1, in2);

  TRLWE c_1 = trlwe_new_sample(in1, key);
  TRLWE_DFT c_dft = trlwe_alloc_new_DFT_sample(k, N);
  TRLWE_DFT res_c = trlwe_alloc_new_DFT_sample(k, N);
  trlwe_to_DFT(c_dft, c_1);
  polynomial_torus_to_DFT(in2_dft, in2);
  trlwe_DFT_mul_by_polynomial(res_c, c_dft, in2_dft);
  trlwe_from_DFT(c_1, res_c);
  trlwe_phase(res, c_1, key);
  TEST_ASSERT_TORUS_ARRAY_WITHIN_MESSAGE(1UL << 58, res_e->coeffs, res->coeffs, N, "TRLWE poly mul failed");

  free_trlwe(c_1);
  free_trlwe(c_dft);
  free_trlwe(res_c);
  free_polynomial(in1);
  free_polynomial(in2);
  free_polynomial(res);
  free_polynomial(res_e);
  free_polynomial(in2_dft);
  free_trlwe_key(key);
}


void test_trlwe_ks(){
  TRLWE_Key key_1 = trlwe_new_binary_key(N, k, rlwe_std_dev);
  TRLWE_Key key_2 = trlwe_new_binary_key(N, k, rlwe_std_dev);
  TRLWE_KS_Key ks_key_1_to_2 = trlwe_new_KS_key(key_2, key_1, 4, 12);
  TRLWE_KS_Key ks_key_2_to_1 = trlwe_new_KS_key(key_1, key_2, 4, 12);

  TorusPolynomial in1 = polynomial_new_torus_polynomial(N);
  TorusPolynomial in2 = polynomial_new_torus_polynomial(N);
  generate_random_bytes(sizeof(Torus)*N, (uint8_t *) in1->coeffs);
  generate_random_bytes(sizeof(Torus)*N, (uint8_t *) in2->coeffs);

  TRLWE c_2 = trlwe_new_sample(in2, key_2);
  TRLWE c_2_k1 = trlwe_alloc_new_sample(k, N);
  TorusPolynomial c_2_k1_r = polynomial_new_torus_polynomial(N);

  trlwe_keyswitch(c_2_k1, c_2, ks_key_2_to_1);
  trlwe_phase(c_2_k1_r, c_2_k1, key_1);
  TEST_ASSERT_TORUS_ARRAY_WITHIN_MESSAGE(1UL << 58, in2->coeffs, c_2_k1_r->coeffs, N, "TRLWE ks_key_2_to_1 failed");

  TRLWE c_1 = trlwe_new_sample(in1, key_1);
  TRLWE c_1_k2 = trlwe_alloc_new_sample(k, N);
  TorusPolynomial c_1_k2_r = polynomial_new_torus_polynomial(N);

  trlwe_keyswitch(c_1_k2, c_1, ks_key_1_to_2);
  trlwe_phase(c_1_k2_r, c_1_k2, key_2);
  TEST_ASSERT_TORUS_ARRAY_WITHIN_MESSAGE(1UL << 58, in1->coeffs, c_1_k2_r->coeffs, N, "TRLWE ks_key_1_to_2 failed");

  free_trlwe(c_1);
  free_trlwe(c_2);
  free_trlwe(c_1_k2);
  free_trlwe(c_2_k1);
  free_trlwe_ks_key(ks_key_1_to_2);
  free_trlwe_ks_key(ks_key_2_to_1);
  free_trlwe_key(key_1);
  free_trlwe_key(key_2);
  free_polynomial(in1);
  free_polynomial(in2);
  free_polynomial(c_1_k2_r);
  free_polynomial(c_2_k1_r);
}


void test_trlwe_full_packing_ks(){
  TLWE_Key key_1 = tlwe_new_binary_key(n, lwe_std_dev);
  TRLWE_Key key_2 = trlwe_new_binary_key(N, k, rlwe_std_dev);
  TRLWE_KS_Key ks_key_1_to_2 = trlwe_new_full_packing_KS_key(key_2, key_1, t, base_bit);

  const uint64_t size = N;
  Torus in[N];
  generate_random_bytes(sizeof(Torus)*N, (uint8_t *) in);
  TLWE in_c[size];
  TorusPolynomial exp = polynomial_new_torus_polynomial(N);
  TorusPolynomial res = polynomial_new_torus_polynomial(N);

  for (size_t i = 0; i < size; i++){
    in_c[i] = tlwe_new_sample(in[i], key_1);
    exp->coeffs[i] = in[i];
  }

  TRLWE out = trlwe_alloc_new_sample(k, N);
  trlwe_full_packing_keyswitch(out, in_c, size, ks_key_1_to_2);
  trlwe_phase(res, out, key_2);
  TEST_ASSERT_TORUS_ARRAY_WITHIN_MESSAGE(1UL << 58, exp->coeffs, res->coeffs, N, "TRLWE ks_key_2_to_1 failed");


  free_trlwe(out);
  for (size_t i = 0; i < size; i++){
    free_tlwe(in_c[i]);
  }
  free_trlwe_ks_key(ks_key_1_to_2);
  free_tlwe_key(key_1);
  free_trlwe_key(key_2);
  free_polynomial(exp);
  free_polynomial(res);
}

void test_trlwe_mul(){
  SKIP_IF_TORUS32
  TRLWE_Key key = trlwe_new_binary_key(N, k, rlwe_std_dev);
  TRLWE_KS_Key rlk = trlwe_new_RL_key(key, 2, 20);
  TLWE_Key tlwe_key = tlwe_alloc_key(N, rlwe_std_dev);
  trlwe_extract_tlwe_key(tlwe_key, key);
  _glb_debug_trlwe_key = key;
  _glb_debug_tlwe_key = tlwe_key;


  TRLWE c_1 = trlwe_new_sample(NULL, key);
  TRLWE c_2 = trlwe_new_sample(NULL, key);
  TRLWE c_out = trlwe_alloc_new_sample(k, N);
  TLWE tlwe_out = tlwe_alloc_sample(N);

  const int precision = 4, mod_mask = (1<<precision) - 1;
  int in1, in2, exp_res, res;
  generate_random_bytes(sizeof(int), (uint8_t *) &in1);
  generate_random_bytes(sizeof(int), (uint8_t *) &in2);
  in1 &= mod_mask;
  in2 &= mod_mask;
  exp_res = (in1*in2)&mod_mask;
  printf("%d*%d=%d\n", in1,in2, exp_res);

  c_1->b->coeffs[0] += int2torus(in1, precision);
  c_2->b->coeffs[0] += int2torus(in2, precision);
  trlwe_tensor_prod_FFT(c_out, c_1, c_2, precision, rlk);
  trlwe_extract_tlwe(tlwe_out, c_out, 0);
  res = torus2int(tlwe_phase(tlwe_out, tlwe_key), precision);
  TEST_ASSERT_EQUAL_INT32_MESSAGE(exp_res, res, "Tensor prod failed.");

  free_trlwe(c_1);
  free_trlwe(c_2);
  free_trlwe(c_out);
  free_trlwe_ks_key(rlk);
  free_trlwe_key(key);
  free_tlwe(tlwe_out);
  free_tlwe_key(tlwe_key);
}

void test_tlwe_mul(){
  SKIP_IF_TORUS32
  TRLWE_Key key = trlwe_new_binary_key(N, k, rlwe_std_dev);
  TRLWE_KS_Key rlk = trlwe_new_RL_key(key, 2, 20);
  TLWE_Key tlwe_key = tlwe_alloc_key(N, rlwe_std_dev);
  trlwe_extract_tlwe_key(tlwe_key, key);
  _glb_debug_trlwe_key = key;
  _glb_debug_tlwe_key = tlwe_key;
  Generic_KS_Key ksk = trlwe_new_packing1_KS_key(key, tlwe_key, t, base_bit);

  const int precision = 4, mod_mask = (1<<precision) - 1;
  int in1, in2, exp_res, res;
  generate_random_bytes(sizeof(int), (uint8_t *) &in1);
  generate_random_bytes(sizeof(int), (uint8_t *) &in2);
  in1 &= mod_mask;
  in2 &= mod_mask;
  TLWE c_1 = tlwe_new_sample(int2torus(in1, precision), tlwe_key);
  TLWE c_2 = tlwe_new_sample(int2torus(in2, precision), tlwe_key);
  TLWE c_out = tlwe_alloc_sample(N);

  tlwe_mul(c_out, c_1, c_2, precision, ksk, rlk);
  
  exp_res = (in1*in2)&mod_mask;
  printf("%d*%d=%d\n", in1,in2, exp_res);

  res = torus2int(tlwe_phase(c_out, tlwe_key), precision);
  TEST_ASSERT_EQUAL_INT32_MESSAGE(exp_res, res, "TLWE mul failed.");

  free_tlwe(c_1);
  free_tlwe(c_2);
  free_tlwe(c_out);
  free_trlwe_ks_key(rlk);
  free_trlwe_key(key);
  free_tlwe_key(tlwe_key);
}


void test_trlwe_packing_ks(){
  TLWE_Key key_tlwe = tlwe_new_binary_key(N, rlwe_std_dev);
  TRLWE_Key key_trlwe = trlwe_new_binary_key(N, k, rlwe_std_dev);
  trlwe_extract_tlwe_key(key_tlwe, key_trlwe);

  
  LUT_Packing_KS_Key ks_key = trlwe_new_packing_KS_key(key_trlwe, key_tlwe, t, base_bit, 4);

  Torus in[4];
  TLWE c[4];
  generate_random_bytes(sizeof(Torus)*4, (uint8_t *) &in);
  for (size_t i = 0; i < 4; i++){
    c[i] = tlwe_new_sample(in[i], key_tlwe);
  }

  TRLWE c_out = trlwe_new_noiseless_trivial_sample(0, k, N);
  trlwe_packing_keyswitch(c_out, c, ks_key);

  TorusPolynomial poly_out = polynomial_new_torus_polynomial(N);
  trlwe_phase(poly_out, c_out, key_trlwe);

  TorusPolynomial poly_res = polynomial_new_torus_polynomial(N);
  for (size_t i = 0; i < N; i++) poly_res->coeffs[i] = in[i/(N/4)];

  TEST_ASSERT_TORUS_ARRAY_WITHIN_MESSAGE(1UL << 58, poly_res->coeffs, poly_out->coeffs, N, "TRLWE Packing KS failed.");

  free_tlwe_key(key_tlwe);
  free_trlwe_key(key_trlwe);
  free_trlwe_packing_ks_key(ks_key);
  for (size_t i = 0; i < 4; i++) free_tlwe(c[i]);
  free_trlwe(c_out);
  free_polynomial(poly_out);
  free_polynomial(poly_res);
}

void test_blind_rotate(){
  TRLWE_Key key_trlwe = trlwe_new_binary_key(N, k, rlwe_std_dev);
  TRGSW_Key trgsw_key = trgsw_new_key(key_trlwe, l, Bg_bit);
  _glb_debug_trgsw_key = trgsw_key;
  _glb_debug_trlwe_key = key_trlwe;

  TRLWE lut = trlwe_new_noiseless_trivial_sample(0, k, N);
  for (size_t i = 0; i < N; i++) lut->b->coeffs[i] = double2torus(1./8);
  lut->b->coeffs[0] = double2torus(2./8); 

  Torus in;
  generate_random_bytes(sizeof(Torus), (uint8_t *) &in);

  TRGSW tmp = trgsw_new_monomial_sample(1, 0, trgsw_key);
  TRGSW_DFT s = trgsw_alloc_new_DFT_sample(l, Bg_bit, k, N);
  trgsw_to_DFT(s, tmp);
  blind_rotate(lut, &in, &s, 1);

  TorusPolynomial poly_out = polynomial_new_torus_polynomial(N);
  trlwe_phase(poly_out, lut, key_trlwe);

  TorusPolynomial poly_res = polynomial_new_torus_polynomial(N);
  TorusPolynomial poly_res_2 = polynomial_new_torus_polynomial(N);
  for (size_t i = 0; i < N; i++) poly_res->coeffs[i] = double2torus(1./8);
  poly_res->coeffs[0] = double2torus(2./8); 

  const int ai = torus2int(in, (int) log2(2*N));
  torus_polynomial_mul_by_xai(poly_res_2, poly_res, ai);

  TEST_ASSERT_TORUS_ARRAY_WITHIN_MESSAGE(1UL << 58, poly_res_2->coeffs, poly_out->coeffs, N, "Blind Rotate failed.");
  free_trgsw(tmp);
  free_trgsw(s);
  free_trlwe(lut);
  free_polynomial(poly_out);
  free_polynomial(poly_res);
  free_polynomial(poly_res_2);
  free_trlwe_key(key_trlwe);
  free_trgsw_key(trgsw_key);
}

void test_functional_bootstrap_unfolded(){
  TLWE_Key key_tlwe = tlwe_new_binary_key(632, lwe_std_dev); // n must be multiple of unfolding
  TLWE_Key key_tlwe_out = tlwe_new_binary_key(N, rlwe_std_dev);
  TRLWE_Key key_trlwe = trlwe_new_binary_key(N, k, rlwe_std_dev);
  trlwe_extract_tlwe_key(key_tlwe_out, key_trlwe);

  TRGSW_Key trgsw_key = trgsw_new_key(key_trlwe, l, Bg_bit);

  Torus in[5];
  TLWE c[6];
  generate_random_bytes(sizeof(Torus)*4, (uint8_t *) in);
  for (size_t i = 0; i < 5; i++){
    c[i] = tlwe_new_sample(in[i], key_tlwe_out);
  }
  c[5] = tlwe_new_sample(double2torus(1./8), key_tlwe);

  TorusPolynomial poly_res = polynomial_new_torus_polynomial(N);
  for (size_t i = 0; i < N; i++) poly_res->coeffs[i] = in[i/(N/4)];

  TRLWE lut_c = trlwe_new_noiseless_trivial_sample(poly_res, k, N);

  Bootstrap_Key bk_key = new_bootstrap_key(trgsw_key, key_tlwe, 2);
  functional_bootstrap(c[4], lut_c, c[5], bk_key, 4);

  TEST_ASSERT_TORUS_WITHIN_MESSAGE(1UL << 58, in[1], tlwe_phase(c[4], key_tlwe_out), "Bootstrap failed. (cleartext LUT, unfolding = 2)");
  
  free_bootstrap_key(bk_key);
  bk_key = new_bootstrap_key(trgsw_key, key_tlwe, 4);
  functional_bootstrap(c[4], lut_c, c[5], bk_key, 4);

  TEST_ASSERT_TORUS_WITHIN_MESSAGE(1UL << 58, in[1], tlwe_phase(c[4], key_tlwe_out), "Bootstrap failed. (cleartext LUT, unfolding = 4)");

  free_bootstrap_key(bk_key);
  bk_key = new_bootstrap_key(trgsw_key, key_tlwe, 8);
  functional_bootstrap(c[4], lut_c, c[5], bk_key, 4);

  // Unfoldings greater than 4 might fail eventually with the predefined (TFHEpp's) parameters
  TEST_ASSERT_TORUS_WITHIN_MESSAGE(1UL << 58, in[1], tlwe_phase(c[4], key_tlwe_out), "Bootstrap failed. (cleartext LUT, unfolding = 8)");

  free_tlwe_key(key_tlwe);
  free_tlwe_key(key_tlwe_out);
  free_trlwe_key(key_trlwe);
  free_trgsw_key(trgsw_key);
  free_bootstrap_key(bk_key);
  for (size_t i = 0; i < 6; i++) free_tlwe(c[i]);
  free_polynomial(poly_res);
  free_trlwe(lut_c);
}


void test_programmable_bootstrap(){
  TLWE_Key key_tlwe = tlwe_new_binary_key(n, lwe_std_dev);
  TLWE_Key key_tlwe_out = tlwe_new_binary_key(N, rlwe_std_dev);
  TRLWE_Key key_trlwe = trlwe_new_binary_key(N, k, rlwe_std_dev);
  trlwe_extract_tlwe_key(key_tlwe_out, key_trlwe);

  TRGSW_Key trgsw_key = trgsw_new_key(key_trlwe, l, Bg_bit);
  Bootstrap_Key bk_key = new_bootstrap_key(trgsw_key, key_tlwe, 1);

  Torus in[5];
  TLWE c[6];
  generate_random_bytes(sizeof(Torus)*4, (uint8_t *) in);
  for (size_t i = 0; i < 5; i++){
    c[i] = tlwe_new_sample(in[i], key_tlwe_out);
  }
  c[5] = tlwe_new_sample(double2torus(1./8), key_tlwe);

  TorusPolynomial poly_res = polynomial_new_torus_polynomial(N);
  for (size_t i = 0; i < N; i++) poly_res->coeffs[i] = in[i/(N/4)];

  TRLWE lut_c = trlwe_new_noiseless_trivial_sample(poly_res, k, N);

  programmable_bootstrap(c[4], lut_c, c[5], bk_key, 3, 0, 0);

  TEST_ASSERT_TORUS_WITHIN_MESSAGE(1UL << 58, in[1], tlwe_phase(c[4], key_tlwe_out), "Bootstrap with cleartext LUT failed 0,0.");

  tlwe_sample(c[5], int2torus(0xA, 6), key_tlwe);
  programmable_bootstrap(c[4], lut_c, c[5], bk_key, 3, 3, 0);

  TEST_ASSERT_TORUS_WITHIN_MESSAGE(1UL << 58, in[2], tlwe_phase(c[4], key_tlwe_out), "Bootstrap with cleartext LUT failed 3,0.");

  free_tlwe_key(key_tlwe);
  free_tlwe_key(key_tlwe_out);
  free_trlwe_key(key_trlwe);
  free_trgsw_key(trgsw_key);
  free_bootstrap_key(bk_key);
  for (size_t i = 0; i < 6; i++) free_tlwe(c[i]);
  free_polynomial(poly_res);
  free_trlwe(lut_c);
}


void test_functional_bootstrap(){
  TLWE_Key key_tlwe = tlwe_new_binary_key(n, lwe_std_dev);
  TLWE_Key key_tlwe_out = tlwe_new_binary_key(N, rlwe_std_dev);
  TRLWE_Key key_trlwe = trlwe_new_binary_key(N, k, rlwe_std_dev);
  trlwe_extract_tlwe_key(key_tlwe_out, key_trlwe);

  TRGSW_Key trgsw_key = trgsw_new_key(key_trlwe, l, Bg_bit);
  Bootstrap_Key bk_key = new_bootstrap_key(trgsw_key, key_tlwe, 1);

  Torus in[5];
  TLWE c[6];
  generate_random_bytes(sizeof(Torus)*4, (uint8_t *) in);
  for (size_t i = 0; i < 5; i++){
    c[i] = tlwe_new_sample(in[i], key_tlwe_out);
  }
  c[5] = tlwe_new_sample(double2torus(1./8), key_tlwe);

  TorusPolynomial poly_res = polynomial_new_torus_polynomial(N);
  for (size_t i = 0; i < N; i++) poly_res->coeffs[i] = in[i/(N/4)];

  TRLWE lut_c = trlwe_new_noiseless_trivial_sample(poly_res, k, N);

  functional_bootstrap(c[4], lut_c, c[5], bk_key, 4);

  TEST_ASSERT_TORUS_WITHIN_MESSAGE(1UL << 58, in[1], tlwe_phase(c[4], key_tlwe_out), "Bootstrap with cleartext LUT failed.");

  free_tlwe_key(key_tlwe);
  free_tlwe_key(key_tlwe_out);
  free_trlwe_key(key_trlwe);
  free_trgsw_key(trgsw_key);
  free_bootstrap_key(bk_key);
  for (size_t i = 0; i < 6; i++) free_tlwe(c[i]);
  free_polynomial(poly_res);
  free_trlwe(lut_c);
}


void test_functional_bootstrap_ga(){
  TLWE_Key key_tlwe = tlwe_new_binary_key(n, lwe_std_dev);
  TLWE_Key key_tlwe_out = tlwe_new_binary_key(N, rlwe_std_dev);
  TRLWE_Key key_trlwe = trlwe_new_binary_key(N, k, rlwe_std_dev);
  trlwe_extract_tlwe_key(key_tlwe_out, key_trlwe);

  TRGSW_Key trgsw_key = trgsw_new_key(key_trlwe, l, Bg_bit);
  Bootstrap_GA_Key bk_key = new_bootstrap_key_ga(trgsw_key, key_tlwe);

  Torus in[5];
  TLWE c[6];
  generate_random_bytes(sizeof(Torus)*4, (uint8_t *) in);
  for (size_t i = 0; i < 5; i++){
    c[i] = tlwe_new_sample(in[i], key_tlwe_out);
  }
  c[5] = tlwe_new_sample(double2torus(1./8), key_tlwe);

  TorusPolynomial poly_res = polynomial_new_torus_polynomial(N);
  for (size_t i = 0; i < N; i++) poly_res->coeffs[i] = in[i/(N/4)];

  TRLWE lut_c = trlwe_new_noiseless_trivial_sample(poly_res, k, N);

  functional_bootstrap_ga(c[4], lut_c, c[5], bk_key, 4);

  TEST_ASSERT_TORUS_WITHIN_MESSAGE(1UL << 58, in[1], tlwe_phase(c[4], key_tlwe_out), "Bootstrap GA with cleartext LUT failed.");

  free_tlwe_key(key_tlwe);
  free_tlwe_key(key_tlwe_out);
  free_trlwe_key(key_trlwe);
  free_trgsw_key(trgsw_key);
  free_bootstrap_key_ga(bk_key);
  for (size_t i = 0; i < 6; i++) free_tlwe(c[i]);
  free_polynomial(poly_res);
  free_trlwe(lut_c);
}


void test_functional_bootstrap_ga_bounded_key(){
  TLWE_Key key_tlwe = tlwe_new_bounded_key(n, 4, lwe_std_dev);
  TLWE_Key key_tlwe_out = tlwe_new_binary_key(N, rlwe_std_dev);
  TRLWE_Key key_trlwe = trlwe_new_bounded_key(N, k, 16, rlwe_std_dev);
  trlwe_extract_tlwe_key(key_tlwe_out, key_trlwe);

  TRGSW_Key trgsw_key = trgsw_new_key(key_trlwe, l, Bg_bit);
  Bootstrap_GA_Key bk_key = new_bootstrap_key_ga(trgsw_key, key_tlwe);

  Torus in[5];
  TLWE c[6];
  generate_random_bytes(sizeof(Torus)*4, (uint8_t *) in);
  for (size_t i = 0; i < 5; i++){
    c[i] = tlwe_new_sample(in[i], key_tlwe_out);
  }
  c[5] = tlwe_new_sample(double2torus(1./8), key_tlwe);

  TorusPolynomial poly_res = polynomial_new_torus_polynomial(N);
  for (size_t i = 0; i < N; i++) poly_res->coeffs[i] = in[i/(N/4)];


  TRLWE lut_c = trlwe_new_noiseless_trivial_sample(poly_res, k, N);

  for (size_t i = 0; i < 4; i++)
  {
    tlwe_sample(c[5], double2torus(((double)i)/8.), key_tlwe);
    functional_bootstrap_ga(c[4], lut_c, c[5], bk_key, 4);
    TEST_ASSERT_TORUS_WITHIN_MESSAGE(1UL << 58, in[i], tlwe_phase(c[4], key_tlwe_out), "Bootstrap GA with cleartext LUT failed.");
  }

  free_tlwe_key(key_tlwe);
  free_tlwe_key(key_tlwe_out);
  free_trlwe_key(key_trlwe);
  free_trgsw_key(trgsw_key);
  free_bootstrap_key_ga(bk_key);
  for (size_t i = 0; i < 6; i++) free_tlwe(c[i]);
  free_polynomial(poly_res);
  free_trlwe(lut_c);
}

void test_functional_bootstrap_with_encrypted_LUT(){
  TLWE_Key key_tlwe = tlwe_new_binary_key(n, lwe_std_dev);
  TLWE_Key key_tlwe_out = tlwe_new_binary_key(N, rlwe_std_dev);
  TRLWE_Key key_trlwe = trlwe_new_binary_key(N, k, rlwe_std_dev);
  trlwe_extract_tlwe_key(key_tlwe_out, key_trlwe);

  TRGSW_Key trgsw_key = trgsw_new_key(key_trlwe, l, Bg_bit);
  Bootstrap_Key bk_key = new_bootstrap_key(trgsw_key, key_tlwe, 1);

  Torus in[5];
  TLWE c[6];
  generate_random_bytes(sizeof(Torus)*4, (uint8_t *) in);
  for (size_t i = 0; i < 5; i++){
    c[i] = tlwe_new_sample(in[i], key_tlwe_out);
  }
  c[5] = tlwe_new_sample(double2torus(1./8), key_tlwe);

  LUT_Packing_KS_Key ks_key = trlwe_new_packing_KS_key(key_trlwe, key_tlwe_out, t, base_bit, 4);
  TRLWE lut = trlwe_new_noiseless_trivial_sample(0, k, N);
  trlwe_packing_keyswitch(lut, c, ks_key);

  TorusPolynomial poly_res = polynomial_new_torus_polynomial(N);
  for (size_t i = 0; i < N; i++) poly_res->coeffs[i] = in[i/(N/4)];

  TorusPolynomial poly_out = polynomial_new_torus_polynomial(N);
  trlwe_phase(poly_out, lut, key_trlwe);

  TEST_ASSERT_TORUS_ARRAY_WITHIN_MESSAGE(1UL << 58, poly_res->coeffs, poly_out->coeffs, N, "TRLWE Packing KS failed.");

  functional_bootstrap(c[4], lut, c[5], bk_key, 4);

  TEST_ASSERT_TORUS_WITHIN_MESSAGE(1UL << 58, in[1], tlwe_phase(c[4], key_tlwe_out), "Bootstrap with encrypted LUT failed.");
  free_tlwe_key(key_tlwe);
  free_tlwe_key(key_tlwe_out);
  free_trlwe_key(key_trlwe);
  free_trgsw_key(trgsw_key);
  free_bootstrap_key(bk_key);
  for (size_t i = 0; i < 6; i++) free_tlwe(c[i]);
  free_polynomial(poly_out);
  free_polynomial(poly_res);
  free_trlwe(lut);
  free_trlwe_packing_ks_key(ks_key);
}



void test_functional_bootstrap_trgsw(){
  TLWE_Key key_tlwe = tlwe_new_binary_key(n, lwe_std_dev);
  TLWE_Key key_tlwe_out = tlwe_new_binary_key(N, rlwe_std_dev);
  TRLWE_Key key_trlwe = trlwe_new_binary_key(N, k, rlwe_std_dev);
  trlwe_extract_tlwe_key(key_tlwe_out, key_trlwe);

  TRGSW_Key trgsw_key = trgsw_new_key(key_trlwe, l, Bg_bit);
  Bootstrap_Key bk_key = new_bootstrap_key(trgsw_key, key_tlwe, 1);

  Torus in[5];
  TLWE c[6];
  generate_random_bytes(sizeof(Torus)*4, (uint8_t *) in);
  for (size_t i = 0; i < 5; i++){
    c[i] = tlwe_new_sample(in[i], key_tlwe_out);
  }
  c[5] = tlwe_new_sample(double2torus(1./8), key_tlwe);

  TorusPolynomial poly_res = polynomial_new_torus_polynomial(N);
  for (size_t i = 0; i < N; i++) poly_res->coeffs[i] = in[i/(N/4)];

  TRLWE lut_c = trlwe_new_noiseless_trivial_sample(poly_res, k, N);

  TRGSW_DFT c_trgsw = trgsw_alloc_new_DFT_sample(l, Bg_bit, k, N);
  functional_bootstrap_trgsw_phase1(c_trgsw, c[5], bk_key, 4);
  functional_bootstrap_trgsw_phase2(c[4], c_trgsw, lut_c);

  TEST_ASSERT_TORUS_WITHIN_MESSAGE(1UL << 60, in[1], tlwe_phase(c[4], key_tlwe_out), "Bootstrap TRGSW with cleartext LUT failed.");

  LUT_Packing_KS_Key ks_key = trlwe_new_packing_KS_key(key_trlwe, key_tlwe_out, t, base_bit, 4);
  TRLWE lut = trlwe_new_noiseless_trivial_sample(0, k, N);
  trlwe_packing_keyswitch(lut, c, ks_key);

  TorusPolynomial poly_out = polynomial_new_torus_polynomial(N);
  trlwe_phase(poly_out, lut, key_trlwe);

  TEST_ASSERT_TORUS_ARRAY_WITHIN_MESSAGE(1UL << 60, poly_res->coeffs, poly_out->coeffs, N, "TRLWE Packing KS failed.");

  functional_bootstrap_trgsw_phase1(c_trgsw, c[5], bk_key, 4);
  functional_bootstrap_trgsw_phase2(c[4], c_trgsw, lut);

  TEST_ASSERT_TORUS_WITHIN_MESSAGE(1UL << 60, in[1], tlwe_phase(c[4], key_tlwe_out), "Bootstrap TRGSW with encrypted LUT failed.");
  free_tlwe_key(key_tlwe);
  free_tlwe_key(key_tlwe_out);
  free_trlwe_key(key_trlwe);
  free_trgsw_key(trgsw_key);
  free_bootstrap_key(bk_key);
  for (size_t i = 0; i < 6; i++) free_tlwe(c[i]);
  free_polynomial(poly_out);
  free_polynomial(poly_res);
  free_trlwe(lut);
  free_trlwe(lut_c);
  free_trlwe_packing_ks_key(ks_key);
  free_trgsw(c_trgsw);
}

void test_functional_mv_bootstrap(){

  TLWE_Key key_tlwe = tlwe_new_binary_key(n, lwe_std_dev);
  TLWE_Key key_tlwe_out = tlwe_new_binary_key(N, rlwe_std_dev);
  TRLWE_Key key_trlwe = trlwe_new_binary_key(N, k, rlwe_std_dev);
  trlwe_extract_tlwe_key(key_tlwe_out, key_trlwe);
  _glb_debug_trlwe_key = key_trlwe;
  _glb_debug_tlwe_key = key_tlwe_out;

  TRGSW_Key trgsw_key = trgsw_new_key(key_trlwe, l, Bg_bit);
  Bootstrap_Key bk_key = new_bootstrap_key(trgsw_key, key_tlwe, 1);

  Torus in[5];
  TLWE c[6];
  generate_random_bytes(sizeof(Torus)*4, (uint8_t *) in);
  for (size_t i = 0; i < 5; i++){
    c[i] = tlwe_new_sample(in[i], key_tlwe_out);
  }
  c[5] = tlwe_new_sample(double2torus(1./8), key_tlwe);

  int lut_int[4] = {1,2,3,0};

  TRLWE * bk_out = trlwe_alloc_new_sample_array(5, k, N);
  multivalue_bootstrap_phase1(bk_out, c[5], bk_key, 4);
  multivalue_bootstrap_phase2(c[4], lut_int, bk_out, 4, 2);

  TEST_ASSERT_TORUS_WITHIN_MESSAGE(1UL << 58, double2torus(((double) lut_int[1])/8), tlwe_phase(c[4], key_tlwe_out), "Multi value Bootstrap with cleartext LUT failed.");
  free_tlwe_key(key_tlwe);
  free_tlwe_key(key_tlwe_out);
  free_trlwe_key(key_trlwe);
  free_trgsw_key(trgsw_key);
  free_bootstrap_key(bk_key);
  for (size_t i = 0; i < 6; i++) free_tlwe(c[i]);
  free_trlwe_array(bk_out, 5);
}

void test_io_priv(){
  TLWE_Key key_tlwe = tlwe_new_binary_key(n, lwe_std_dev);
  TRLWE_Key key_trlwe = trlwe_new_binary_key(N, k, rlwe_std_dev);
  TLWE_Key extracted_key = tlwe_alloc_key(N*k, rlwe_std_dev);
  trlwe_extract_tlwe_key(extracted_key, key_trlwe);
  TRGSW_Key trgsw_key = trgsw_new_key(key_trlwe, l, Bg_bit);

  FILE * fd = fopen("priv.key", "w+");
  tlwe_save_key(fd, key_tlwe);
  tlwe_save_key(fd, extracted_key);
  trlwe_save_key(fd, key_trlwe);
  trgsw_save_key(fd, trgsw_key);
  fclose(fd);

  fd = fopen("priv.key", "r");
  TLWE_Key loaded_tlwe_key = tlwe_load_new_key(fd);
  TEST_ASSERT_EQUAL_INT(n, loaded_tlwe_key->n);
  TEST_ASSERT_EQUAL_DOUBLE(lwe_std_dev, loaded_tlwe_key->sigma);
  TEST_ASSERT_EQUAL_CHAR_ARRAY(key_tlwe->s, loaded_tlwe_key->s, n);

  TLWE_Key loaded_extracted_key = tlwe_load_new_key(fd);
  TEST_ASSERT_EQUAL_INT(N*k, loaded_extracted_key->n);
  TEST_ASSERT_EQUAL_DOUBLE(rlwe_std_dev, loaded_extracted_key->sigma);
  TEST_ASSERT_EQUAL_CHAR_ARRAY(extracted_key->s, loaded_extracted_key->s, N*k);

  TRLWE_Key loaded_trlwe_key = trlwe_load_new_key(fd);
  TEST_ASSERT_EQUAL_INT(k, loaded_trlwe_key->k);
  TEST_ASSERT_EQUAL_INT(N, loaded_trlwe_key->s[0]->N);
  TEST_ASSERT_EQUAL_DOUBLE(rlwe_std_dev, loaded_trlwe_key->sigma);
  for (size_t i = 0; i < k; i++){
    TEST_ASSERT_EQUAL_CHAR_ARRAY(key_trlwe->s[i]->coeffs, loaded_trlwe_key->s[i]->coeffs, N);
  }
  

  TRGSW_Key loaded_trgsw_key = trgsw_load_new_key(fd);
  TEST_ASSERT_EQUAL_INT(Bg_bit, loaded_trgsw_key->Bg_bit);
  TEST_ASSERT_EQUAL_INT(l, loaded_trgsw_key->l);
  TEST_ASSERT_EQUAL_INT(k, loaded_trgsw_key->trlwe_key->k);
  TEST_ASSERT_EQUAL_INT(N, loaded_trgsw_key->trlwe_key->s[0]->N);
  TEST_ASSERT_EQUAL_DOUBLE(rlwe_std_dev, loaded_trgsw_key->trlwe_key->sigma);
  for (size_t i = 0; i < k; i++){
    TEST_ASSERT_EQUAL_CHAR_ARRAY(key_trlwe->s[i]->coeffs, loaded_trgsw_key->trlwe_key->s[i]->coeffs, N);
  }

  fclose(fd);
  remove("priv.key");
  free_tlwe_key(key_tlwe);
  free_tlwe_key(extracted_key);
  free_trlwe_key(key_trlwe);
  free_trgsw_key(trgsw_key);
  free_tlwe_key(loaded_tlwe_key);
  free_tlwe_key(loaded_extracted_key);
  free_trlwe_key(loaded_trlwe_key);
  free_trgsw_key(loaded_trgsw_key);
}

void test_io_pub(){
  TLWE_Key key_tlwe = tlwe_new_binary_key(n, lwe_std_dev);
  TRLWE_Key key_trlwe = trlwe_new_binary_key(N, k, rlwe_std_dev);
  TLWE_Key extracted_key = tlwe_alloc_key(N*k, rlwe_std_dev);
  trlwe_extract_tlwe_key(extracted_key, key_trlwe);
  TRGSW_Key trgsw_key = trgsw_new_key(key_trlwe, l, Bg_bit);

  Bootstrap_Key bk_key = new_bootstrap_key(trgsw_key, key_tlwe, 1);
  LUT_Packing_KS_Key packing_ks_key = trlwe_new_packing_KS_key(key_trlwe, extracted_key, t, base_bit, 4);
  TLWE_KS_Key ks_key = tlwe_new_KS_key(key_tlwe, extracted_key, t, base_bit);

  FILE * fd = fopen("pub.key", "w+");
  save_bootstrap_key(fd, bk_key);
  tlwe_save_KS_key(fd, ks_key);
  trlwe_save_packing_KS_key(fd, packing_ks_key);
  fclose(fd);

  // free_trlwe_packing_ks_key(packing_ks_key); 

  fd = fopen("pub.key", "r");

  Bootstrap_Key loaded_bk_key = load_new_bootstrap_key(fd);
  TEST_ASSERT_EQUAL_INT(bk_key->n, loaded_bk_key->n);
  TEST_ASSERT_EQUAL_INT(bk_key->N, loaded_bk_key->N);
  TEST_ASSERT_EQUAL_INT(bk_key->l, loaded_bk_key->l);
  TEST_ASSERT_EQUAL_INT(bk_key->k, loaded_bk_key->k);
  TEST_ASSERT_EQUAL_INT(bk_key->Bg_bit, loaded_bk_key->Bg_bit);
  for (size_t i = 0; i < bk_key->n; i++){
    for (size_t j = 0; j < bk_key->l * (bk_key->k + 1); j++){
      for (size_t q = 0; q < bk_key->k; q++){
        TEST_ASSERT_EQUAL_DOUBLE_ARRAY(bk_key->s[i]->samples[j]->a[q]->coeffs, loaded_bk_key->s[i]->samples[j]->a[q]->coeffs, bk_key->N);
      }      
      TEST_ASSERT_EQUAL_DOUBLE_ARRAY(bk_key->s[i]->samples[j]->b->coeffs, loaded_bk_key->s[i]->samples[j]->b->coeffs, bk_key->N);
    }
  }

  TLWE_KS_Key loaded_ks_key = tlwe_load_new_KS_key(fd);
  TEST_ASSERT_EQUAL_INT(ks_key->n, loaded_ks_key->n);
  TEST_ASSERT_EQUAL_INT(ks_key->t, loaded_ks_key->t);
  TEST_ASSERT_EQUAL_INT(ks_key->base_bit, loaded_ks_key->base_bit);

  for (size_t i = 0; i < ks_key->n; i++){
    for (size_t j = 0; j < ks_key->t; j++){
      for (size_t k = 0; k < (1 << ks_key->base_bit) - 1; k++){
        TEST_ASSERT_EQUAL_HEX64_ARRAY(ks_key->s[i][j][k]->a, loaded_ks_key->s[i][j][k]->a, ks_key->s[i][j][k]->n);
        TEST_ASSERT_EQUAL_HEX64(ks_key->s[i][j][k]->b, loaded_ks_key->s[i][j][k]->b);
      }
    }
  }

  LUT_Packing_KS_Key loaded_packing_ks_key = trlwe_load_new_packing_KS_key(fd);
  TEST_ASSERT_EQUAL_INT(packing_ks_key->base_bit, loaded_packing_ks_key->base_bit);
  TEST_ASSERT_EQUAL_INT(packing_ks_key->n, loaded_packing_ks_key->n);
  TEST_ASSERT_EQUAL_INT(packing_ks_key->t, loaded_packing_ks_key->t);
  TEST_ASSERT_EQUAL_INT(packing_ks_key->torus_base, loaded_packing_ks_key->torus_base);

  for (size_t i = 0; i < packing_ks_key->n; i++){
    for (size_t e = 0; e < packing_ks_key->torus_base; e++){
      for (size_t j = 0; j < packing_ks_key->t; j++){
        for (size_t k = 0; k <  (1 << packing_ks_key->base_bit) - 1; k++){
          for (size_t l = 0; l < packing_ks_key->s[0][0][0][0]->k; l++){
            TEST_ASSERT_EQUAL_HEX64_ARRAY(packing_ks_key->s[i][e][j][k]->a[l]->coeffs, loaded_packing_ks_key->s[i][e][j][k]->a[l]->coeffs, packing_ks_key->s[0][0][0][0]->b->N);
          }
          TEST_ASSERT_EQUAL_HEX64_ARRAY(packing_ks_key->s[i][e][j][k]->b->coeffs, loaded_packing_ks_key->s[i][e][j][k]->b->coeffs, packing_ks_key->s[0][0][0][0]->b->N);
        }
      }
    }
  }

  fclose(fd);
  remove("pub.key");

  free_tlwe_key(key_tlwe);
  free_tlwe_key(extracted_key);
  free_trlwe_key(key_trlwe);
  free_trgsw_key(trgsw_key);
  free_bootstrap_key(bk_key);
  free_tlwe_ks_key(ks_key);
  free_trlwe_packing_ks_key(packing_ks_key); 
  free_bootstrap_key(loaded_bk_key);
  free_tlwe_ks_key(loaded_ks_key);
  free_trlwe_packing_ks_key(loaded_packing_ks_key); 
}

int main(int argc, char const *argv[])
{
  UNITY_BEGIN();
  // RUN_TEST(test_io_pub);
  // RUN_TEST(test_io_priv);
  // RUN_TEST(test_FDFB_CLOT21_3);
  RUN_TEST(test_normal_generator);
  RUN_TEST(test_tlwe);
  RUN_TEST(test_tlwe_ks);
  RUN_TEST(test_poly_DFT);
  RUN_TEST(test_poly_DFT_mul);
  RUN_TEST(test_tlwe_pack1_ks_CDKS21);
  RUN_TEST(test_trlwe);
  RUN_TEST(test_trgsw);
  RUN_TEST(test_trgsw_sub);
  RUN_TEST(test_trgsw_mul_by_xai);
  RUN_TEST(test_trgsw_dft);
  RUN_TEST(test_trgsw_trlwe_mul);
  RUN_TEST(test_blind_rotate);
  RUN_TEST(test_functional_bootstrap);
  RUN_TEST(test_trgsw_mul);
  RUN_TEST(test_programmable_bootstrap);
  RUN_TEST(test_functional_bootstrap_ga);
  RUN_TEST(test_functional_bootstrap_ga_bounded_key);
  RUN_TEST(test_FDFB_KS21);
  RUN_TEST(test_FDFB_new);
  RUN_TEST(test_FDFB_CLOT21);
  RUN_TEST(test_FDFB_CLOT21_2);
  RUN_TEST(test_tlwe_mul);
  RUN_TEST(test_trlwe_mul); // tensor prod
  RUN_TEST(test_trlwe_poly_mul);
  RUN_TEST(test_trlwe_ks);
  RUN_TEST(test_functional_mv_bootstrap);
  RUN_TEST(test_public_mux);
  RUN_TEST(test_compressed_trlwe);
  RUN_TEST(test_compressed_trlwe_rotate_vaes);
  RUN_TEST(test_trlwe_full_packing_ks);
  RUN_TEST(test_multivalue_bootstrap_CLOT21);
  RUN_TEST(test_trlwe_packing_ks);
  RUN_TEST(test_circuit_bootstrap);
  RUN_TEST(test_functional_bootstrap_with_encrypted_LUT);
  // RUN_TEST(test_poly_int128_mul);
  RUN_TEST(test_tlwe_pack_key_priv_ks);
  RUN_TEST(test_tlwe_pack1_ks);
  RUN_TEST(test_trgsw_reg_sub);
  RUN_TEST(test_functional_bootstrap_trgsw);
  RUN_TEST(test_functional_bootstrap_unfolded);
  RUN_TEST(test_trlwe_pack_key_priv_ks);
  // HRD tests
  // RUN_TEST(test_trgsw_reg_sub3);
  // RUN_TEST(test_trgsw_naive_mul);
  // RUN_TEST(test_trgsw_reg_sub2);
  return UNITY_END();
}
