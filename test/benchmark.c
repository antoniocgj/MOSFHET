#include "mosfhet.h"
#include <sys/time.h>
// #include <gperftools/profiler.h>

uint64_t get_time() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return (tv.tv_usec) + (tv.tv_sec * 1000000);
}


#define MAX_EXECS 10000
uint64_t __g_clock_begin, __g_clock_end, __g_clock_array[MAX_EXECS];

void print_bench(char * msg, int n){
  uint64_t mean = 0, sq_err = 0;
  for (size_t i = 0; i < n; i++) mean += __g_clock_array[i];
  mean /= n;
  for (size_t i = 0; i < n; i++) sq_err += (__g_clock_array[i] - mean)*(__g_clock_array[i] - mean);
  double stddev = sqrt(((double) sq_err)/n);
  printf("%s: |%lu,%03lu,%03lu|μs +- |%lf\n", msg, mean/1000000, (mean/1000)%1000, mean%1000, stddev);
}

#define MEASURE_TIME(NAME, REP, MSG, CODE) \
  for (size_t ___i = 0; ___i < REP; ___i++){\
    __g_clock_begin = get_time(); \
    CODE;\
    __g_clock_end = get_time(); \
    __g_clock_array[___i] = __g_clock_end - __g_clock_begin; \
  }\
  print_bench(MSG, REP);
  
// Parameters
// Note: BR Unfolding requires n to be divisible by the unfolding value.
#ifdef TORUS32
// LWE params
const int n = 632;
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
// From Concrete/eprint 2022/704 table 4
#define SET_3
// set 1
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

// Benchmark params
#define _EXECS 100

// Benchmark functions (some depend on others)
// #define BENCH_PRIV_KS
// #define BENCH_PACK1_KS
// #define BENCH_LUT_KS
#define BENCH_MV_BOOTSTRAP
// #define BENCH_TRGSW_BOOTSTRAP
// #define BENCH_CIRCUIT_BOOTSTRAP
// #define BENCH_TENSOR_PROD
// #define BENCH_FDFB
// #define BENCH_UNFOLDING
// #define BENCH_BOOTSTRAP_GA

int main(int argc, char const *argv[]){
  TLWE_Key key_tlwe = tlwe_new_binary_key(n, lwe_std_dev);
  TLWE_Key key_tlwe_out = tlwe_new_binary_key(N, rlwe_std_dev);
  TRLWE_Key key_trlwe = trlwe_new_binary_key(N, k, rlwe_std_dev);
  trlwe_extract_tlwe_key(key_tlwe_out, key_trlwe);
  TRGSW_Key trgsw_key = trgsw_new_key(key_trlwe, l, Bg_bit);
  TLWE_KS_Key tlwe_ksk = tlwe_new_KS_key(key_tlwe, key_tlwe_out, t, base_bit); 

  Torus in[_EXECS*4 + 1];
  TLWE c[_EXECS*4 + 1];
  generate_random_bytes(sizeof(Torus)*(_EXECS*4 + 1), (uint8_t *) in);
  for (size_t i = 0; i < _EXECS*4 + 1; i++) c[i] = tlwe_new_sample(in[i], key_tlwe_out);

  TLWE sel = tlwe_new_sample(double2torus(1./8), key_tlwe);
  TorusPolynomial poly_res = polynomial_new_torus_polynomial(N);
  for (size_t i = 0; i < N; i++) poly_res->coeffs[i] = in[i/(N/4)];

  TRLWE lut_c = trlwe_new_noiseless_trivial_sample(poly_res, k, N);
  Bootstrap_Key bk_key = new_bootstrap_key(trgsw_key, key_tlwe, 1);
  Bootstrap_GA_Key bk_ga_key = new_bootstrap_key_ga(trgsw_key, key_tlwe);

#ifdef BENCH_TRGSW_BOOTSTRAP
  TRGSW_DFT c_trgsw = trgsw_alloc_new_DFT_sample(l, Bg_bit, k, N);
  MEASURE_TIME("TRGSW_BS_P1", _EXECS, "TRGSW Functional Bootstrap PHASE 1",
    functional_bootstrap_trgsw_phase1(c_trgsw, sel, bk_key, 4);
  );

  MEASURE_TIME("TRGSW_BS_P2", _EXECS, "TRGSW Functional Bootstrap PHASE 2",
    functional_bootstrap_trgsw_phase2(c[4], c_trgsw, lut_c);
  );
#endif

#ifdef BENCH_MV_BOOTSTRAP
  int lut_int[4] = {1,2,3,4};
  TRLWE * bk_out = trlwe_alloc_new_sample_array(5, k, N);

  MEASURE_TIME("MVFB_P1", _EXECS, "MV Functional Bootstrap PHASE 1",
    multivalue_bootstrap_phase1(bk_out, sel, bk_key, 4);
  );
  MEASURE_TIME("MVFB_P2", _EXECS, "MV Functional Bootstrap PHASE 2",
    multivalue_bootstrap_phase2(c[4], lut_int, bk_out, 4, 2);
  );
#endif

  TRLWE lut = trlwe_new_noiseless_trivial_sample(0, k, N);

#ifdef BENCH_LUT_KS
  generate_random_bytes(sizeof(Torus)*N, (uint8_t *) lut->b->coeffs);
  LUT_Packing_KS_Key ks_key;

  MEASURE_TIME("ks_keygen", 1, "TRLWE Packing KS KeyGEN",
    ks_key = trlwe_new_packing_KS_key(key_trlwe, key_tlwe_out, t, base_bit, 4);
  );

  MEASURE_TIME("LUT_KS", _EXECS, "TRLWE Packing KS",
    trlwe_packing_keyswitch(lut, &c[___i*4], ks_key);
  );
#endif
 

#ifdef BENCH_PRIV_KS
  Generic_KS_Key kska = trlwe_new_priv_SK_KS_key(key_trlwe, key_tlwe_out, t, base_bit);
  generate_random_bytes(sizeof(Torus)*N, (uint8_t *) lut->b->coeffs);

  MEASURE_TIME("PRIV_KS", _EXECS, "TRLWE Priv KS",
    trlwe_priv_keyswitch(lut, c[___i], kska);
  );

#endif

#ifdef BENCH_PACK1_KS

  Generic_KS_Key kskb = trlwe_new_packing1_KS_key(key_trlwe, key_tlwe_out, t, base_bit);
  MEASURE_TIME("PACK1_KS", _EXECS, "TRLWE Packing 1 KS",
    trlwe_packing1_keyswitch(lut, c[___i], kskb);
  );

  MEASURE_TIME("TLWE_KS", _EXECS, "TLWE KS",
    tlwe_keyswitch(sel, c[___i], tlwe_ksk);
  );

#endif

#ifdef BENCH_CDKS21_KS

  TRLWE_KS_Key *cdks21_ksk = trlwe_new_packing1_KS_key_CDKS21(key_trlwe, key_tlwe_out, l, Bg_bit);

  MEASURE_TIME("TRLWE_KS_CDKS21", _EXECS, "TRLWE Packing 1 KS CDKS21",
    trlwe_packing1_keyswitch_CDKS21(lut, c[___i], cdks21_ksk);
  );

#endif

#ifdef BENCH_CIRCUIT_BOOTSTRAP
  TRGSW trgsw_out = trgsw_alloc_new_sample(l, Bg_bit, k, N);

  MEASURE_TIME("CB", _EXECS, "Circuit Bootstrap",
    circuit_bootstrap(trgsw_out, sel, bk_key, kska, kskb);
  );

  MEASURE_TIME("CB_OPT", _EXECS, "Circuit Bootstrap using Multi value Bootstrap",
    circuit_bootstrap_2(trgsw_out, sel, bk_key, kska, kskb);
  );

#endif

#ifdef BENCH_TENSOR_PROD

  TRLWE c_1 = trlwe_new_sample(NULL, key_trlwe);
  TRLWE c_2 = trlwe_new_sample(NULL, key_trlwe);
  TRLWE c_3 = trlwe_new_sample(NULL, key_trlwe);
  TRLWE_KS_Key rlk = trlwe_new_RL_key(key_trlwe, 2, 20);
  MEASURE_TIME("TENSOR_PROD", _EXECS, "Tensor Product Karatsuba",
    trlwe_tensor_prod(c_3, c_1, c_2, 4, rlk);
  );

  MEASURE_TIME("TENSOR_PROD_FFT", _EXECS, "Tensor Product FFT",
    trlwe_tensor_prod_FFT(c_3, c_1, c_2, 4, rlk);
  );

#endif

#ifdef BENCH_FDFB
  TorusPolynomial poly_res2 = polynomial_new_torus_polynomial(2*N);
  for (size_t i = 0; i < 2*N; i++) poly_res2->coeffs[i] = in[i/(N/4)];

  MEASURE_TIME("FDFB_KS21_OPT", _EXECS, "Full Domain Functional Bootstrap KS21 OPT",
    full_domain_functional_bootstrap_KS21(c[4], poly_res2, sel, bk_key, kskb, 8);
  );

  MEASURE_TIME("FDFB_KS21", _EXECS, "Full Domain Functional Bootstrap KS21",
    full_domain_functional_bootstrap_KS21_2(c[4], poly_res2, sel, bk_key, kskb, 8);
  );

  TRLWE tv[2] = {trlwe_new_noiseless_trivial_sample(NULL, k, N),
                 trlwe_new_noiseless_trivial_sample(NULL, k, N)};
  generate_random_bytes(sizeof(Torus)*N, (uint8_t *) tv[0]->b->coeffs);
  generate_random_bytes(sizeof(Torus)*N, (uint8_t *) tv[1]->b->coeffs);
  
  MEASURE_TIME("FDFB_CLOT21_OPT", _EXECS, "Full Domain Functional Bootstrap CLOT21 OPT",
    full_domain_functional_bootstrap_CLOT21_2(c[4], in, sel, bk_key, kskb, rlk, 4);
  );

  MEASURE_TIME("FDFB_CLOT21", _EXECS, "Full Domain Functional Bootstrap CLOT21",
    full_domain_functional_bootstrap_CLOT21(c[4], tv, sel, bk_key, kskb, rlk, 4);
  );

  MEASURE_TIME("FDFB", _EXECS, "Full Domain Functional Bootstrap (this work)",
    full_domain_functional_bootstrap(c[4], tv[0], sel, bk_key, tlwe_ksk, 3);
  );

#endif

#ifdef BENCH_BOOTSTRAP_GA
  MEASURE_TIME("BK_GA", _EXECS, "Functional Bootstrap GA",
    functional_bootstrap_ga(c[4], lut, sel, bk_ga_key, 4);
  );
  free_bootstrap_key_ga(bk_ga_key);
#endif


#ifdef BENCH_UNFOLDING
  MEASURE_TIME("FB_U1", _EXECS, "Functional Bootstrap Unfold=1",
    functional_bootstrap(c[4], lut, sel, bk_key, 4);
  );
  free_bootstrap_key(bk_key);
  bk_key = new_bootstrap_key(trgsw_key, key_tlwe, 2);

  MEASURE_TIME("FB_U2", _EXECS, "Functional Bootstrap Unfold=2",
    functional_bootstrap(c[4], lut, sel, bk_key, 4);
  );
  free_bootstrap_key(bk_key);
  bk_key = new_bootstrap_key(trgsw_key, key_tlwe, 4);

  MEASURE_TIME("FB_U4", _EXECS, "Functional Bootstrap Unfold=4",
    functional_bootstrap(c[4], lut, sel, bk_key, 4);
  );

  free_bootstrap_key(bk_key);
  bk_key = new_bootstrap_key(trgsw_key, key_tlwe, 8);
  MEASURE_TIME("FB_U8", _EXECS, "Functional Bootstrap Unfold=8",
    functional_bootstrap(c[4], lut, sel, bk_key, 4);
  );

#endif

  return 0;
}
